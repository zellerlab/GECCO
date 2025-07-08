import argparse
import csv
import collections
import itertools
import os
import operator
import pathlib
import random
import typing
from typing import List, Optional, Tuple, Iterable, Callable, Type

import rich.progress
from rich.console import Console

from . import _parser, _common
from .._log import ConsoleLogger


def configure_parser(
    parser: argparse.ArgumentParser, 
    console: Console,
    program: str,
    version: str,
):
    _parser.configure_common(parser, console, program, version)

    _parser.configure_group_input_tables(parser)
    _parser.configure_group_domain_filter(parser)
    _parser.configure_group_training_data(parser)
    _parser.configure_group_training_parameters(parser)

    params_cv = parser.add_argument_group("Cross-Validation")
    params_cv.add_argument(
        "--loto",
        action="store_true",
        default=False,
        help=(
            "Use Leave-One-Type-Out (LOTO) cross-validation instead of "
            "K-folds cross-validation."
        ),
    )
    params_cv.add_argument(
        "--splits",
        type=int,
        default=10,
        help=("Number of folds for cross-validation (if running K-folds)."),
    )

    params_output = parser.add_argument_group(
        "Output",
    )
    params_output.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("cv.tsv"),
        help=(
            "The name of the output file where the cross-validation table "
            "will be written."
        ),
    )

    parser.set_defaults(run=run)


def _group_genes(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    shuffle: bool,
) -> List[List["Gene"]]:
    logger.info("Grouping", "genes by source sequence")
    groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
    seqs = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]
    if shuffle:
        logger.info("Shuffling", "training data sequences")
        random.shuffle(seqs)
    return seqs


def _loto_splits(
    logger: ConsoleLogger, seqs: List[List["Gene"]], *, clusters: pathlib.Path
) -> List[Tuple["NDArray[numpy.bool_]", "NDArray[numpy.bool_]"]]:
    from ...crf.cv import LeaveOneGroupOut
    from ...model import ClusterTable, ClusterType

    logger.info("Loading", "the clusters table")
    with zopen(clusters) as in_:
        table = ClusterTable.load(in_)
        index = {row.sequence_id: row.type for row in table}
        if len(index) != len(table):
            raise ValueError("Training data contains several clusters per sequence")

    logger.info("Grouping", "sequences by cluster types")
    groups = []
    for cluster in seqs:
        ty = next((index[g.source.id] for g in cluster if g.source.id in index), None)
        if ty is None:
            seq_id = next(gene.source.id for gene in cluster)
            logger.warn("Failed", f"to find type of cluster in {seq_id!r}")
            ty = ClusterType()
        groups.append(ty.unpack())

    return list(LeaveOneGroupOut().split(seqs, groups=groups))  # type: ignore


def _kfold_splits(
    logger: ConsoleLogger,
    seqs: List[List["Gene"]],
    *,
    splits: int,
) -> List[Tuple["NDArray[numpy.bool_]", "NDArray[numpy.bool_]"]]:
    import sklearn.model_selection

    return list(sklearn.model_selection.KFold(splits).split(seqs))


def _get_train_data(
    train_indices: Iterable[int], seqs: List[List["Gene"]]
) -> List["Gene"]:
    # extract train data
    return [gene for i in train_indices for gene in seqs[i]]


def _get_test_data(
    test_indices: Iterable[int], seqs: List[List["Gene"]]
) -> List["Gene"]:
    # make a clean copy of the test data without gene probabilities
    return [
        gene.with_protein(
            gene.protein.with_domains(
                d.with_probability(None) for d in gene.protein.domains
            )
        )
        for i in test_indices
        for gene in seqs[i]
    ]


def _fit_predict(
    logger,
    train_data: List["Gene"],
    test_data: List["Gene"],
    *,
    feature_type: str,
    c1: float,
    c2: float,
    window_size: int,
    window_step: int,
    shuffle: bool,
    select: Optional[float],
    correction: Optional[str],
    jobs: int = 0,
) -> List["Gene"]:
    from ...crf import ClusterCRF

    crf = _common.fit_model(
        logger,
        train_data,
        feature_type=feature_type,
        c1=c1,
        c2=c2,
        window_size=window_size,
        window_step=window_step,
        shuffle=shuffle,
        select=select,
        correction=correction,
        jobs=jobs,
    )
    return crf.predict_probabilities(test_data)


def _write_fold(
    logger,
    fold: int,
    truth: List["Genes"],
    predicted: List["Gene"],
    output: pathlib.Path,
    append: bool = False,
) -> None:
    import polars
    from ...model import GeneTable

    frame = (
        GeneTable.from_genes(predicted)
        .data.with_columns(polars.lit(fold).alias("fold"))
        .with_columns(
            GeneTable.from_genes(truth).data.select(
                is_cluster=polars.col("average_p") > 0.5
            )
        )
    )
    with open(output, "ab" if append else "wb") as out:
        frame.write_csv(out, include_header=not append, separator="\t")


def _report_fold(
    logger: ConsoleLogger,
    fold: Optional[int],
    truth: List["Genes"],
    predicted: List["Genes"],
) -> None:
    from sklearn.metrics import average_precision_score, roc_auc_score

    probas = [gene.average_probability for gene in predicted]
    labels = [gene.average_probability > 0.5 for gene in truth]

    aupr = average_precision_score(labels, probas)
    auroc = roc_auc_score(labels, probas)
    if fold:
        logger.info(
            f"Finished training fold {fold} (AUROC={auroc:.3f}, AUPR={aupr:.3f})"
        )
    else:
        logger.info(f"Finished cross validation (AUROC={auroc:.3f}, AUPR={aupr:.3f})")


def run(
    args: argparse.Namespace,
    logger: ConsoleLogger,
    crf_type: Type["ClusterCRF"],
    classifier_type: Type["TypeClassifier"],
    default_hmms: Callable[[], Iterable["HMM"]],
) -> int:
    # seed RNG
    _common.seed_rng(logger, args.seed)

    # load features
    genes = list(_common.load_genes(logger, args.genes))
    features = _common.load_features(logger, args.features)

    # label genes
    genes = _common.annotate_genes(logger, genes, features)

    # load clusters and label genes inside clusters
    clusters = _common.load_clusters(logger, args.clusters)
    genes = _common.label_genes(logger, genes, clusters)

    # group genes
    seqs = _group_genes(logger, genes, shuffle=args.shuffle)
    logger.success("Grouped", "genes into", len(seqs), "sequences")

    # split CV folds
    if args.loto:
        splits = _loto_splits(logger, seqs, clusters=args.clusters)
    else:
        splits = _kfold_splits(logger, seqs, splits=args.splits)

    # run CV
    with logger.progress() as progress:
        unit = "fold" if len(splits) == 1 else "folds"
        task = progress.add_task(
            description="Cross-validating", total=len(splits), unit=unit, precision=""
        )
        logger.info("Performing cross-validation")
        predicted = []
        for i, (train_indices, test_indices) in enumerate(
            progress.track(splits, task_id=task)
        ):
            train_data = _get_train_data(train_indices, seqs)
            test_data = _get_test_data(test_indices, seqs)
            old_genes = _get_train_data(test_indices, seqs)
            new_genes = _fit_predict(
                logger,
                train_data,
                test_data,
                feature_type=args.feature_type,
                c1=args.c1,
                c2=args.c2,
                window_size=args.window_size,
                window_step=args.window_step,
                shuffle=args.shuffle,
                select=args.select,
                correction=args.correction,
                jobs=args.jobs,
                crf_type=crf_type,
            )
            _write_fold(
                logger, i + 1, old_genes, new_genes, output=args.output, append=i > 0
            )
            _report_fold(logger, i + 1, old_genes, new_genes)
            predicted.extend(new_genes)
        _report_fold(logger, None, genes, predicted)
