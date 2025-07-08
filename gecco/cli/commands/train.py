import argparse
import csv
import collections
import itertools
import os
import operator
import pathlib
from typing import List, Optional, Type, Callable, Iterable

from rich.console import Console

from . import _parser, _common
from .._log import ConsoleLogger


def configure_parser(parser: argparse.ArgumentParser, console: Console):
    parser.add_argument(
        "-h", 
        "--help", 
        action=_parser.ConsoleHelpAction,
        help="Show this help message and exit.",
        console=console,
    )

    params_arguments = _parser.configure_group_input_tables(parser)

    params_parameters = parser.add_argument_group(
        "Parameters",
    )
    params_parameters.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=0,
        help=(
            "The number of jobs to use for multithreading. Use 0 to use all "
            "available CPUs."
        ),
    )

    group = parser.add_argument_group("Domain Annotation")
    group.add_argument(
        "-e",
        "--e-filter",
        type=float,
        default=None,
        help=(
            "The e-value cutoff for protein domains to be included. This is "
            "not stable across versions, so consider using a p-value filter "
            "instead."
        ),
    )
    group.add_argument(
        "-p",
        "--p-filter",
        type=float,
        default=1e-9,
        help=("The p-value cutoff for protein domains to be included."),
    )

    group = _parser.configure_group_training_data(parser)
    group = _parser.configure_group_training_parameters(parser)

    group = parser.add_argument_group(
        "Output",
    )
    group.add_argument(
        "-o",
        "--output-dir",
        help="The directory in which to write the output files.",
        type=pathlib.Path,
        default=".",
    )

    parser.set_defaults(run=run)


def _save_model(
    logger: ConsoleLogger, crf: "ClusterCRF", output_dir: pathlib.Path
) -> None:
    model_out = os.path.join(output_dir, "model.pkl")
    logger.info("Pickling", "the model to", repr(model_out))
    with open(model_out, "wb") as out:
        pickle.dump(crf, out, protocol=4)

    logger.info("Computing", "pickled model checksum", level=2)
    hasher = hashlib.md5()
    with open(model_out, "rb") as out:
        for chunk in iter(lambda: out.read(io.DEFAULT_BUFFER_SIZE), b""):
            hasher.update(chunk)

    logger.info(
        "Writing", "pickled model checksum to", repr(f"{model_out}.md5"), level=2
    )
    with open(f"{model_out}.md5", "w") as out_hash:
        out_hash.write(hasher.hexdigest())


def _save_transitions(
    logger: ConsoleLogger, crf: "ClusterCRF", output_dir: pathlib.Path
) -> None:
    logger.info("Writing", "CRF transitions weights")
    with open(os.path.join(output_dir, "model.trans.tsv"), "w") as f:
        writer = csv.writer(f, dialect="excel-tab")
        writer.writerow(["from", "to", "weight"])
        for labels, weight in crf.model.transition_features_.items():
            writer.writerow([*labels, weight])


def _save_weights(
    logger: ConsoleLogger, crf: "ClusterCRF", output_dir: pathlib.Path
) -> None:
    logger.info("Writing", "state weights")
    with open(os.path.join(output_dir, "model.state.tsv"), "w") as f:
        writer = csv.writer(f, dialect="excel-tab")
        writer.writerow(["attr", "label", "weight"])
        for attrs, weight in crf.model.state_features_.items():
            writer.writerow([*attrs, weight])


def _assign_clusters(
    logger: ConsoleLogger, genes: List["Gene"], clusters: "ClusterTable"
) -> List["Cluster"]:
    from ...model import Cluster, ClusterType

    cluster_types = {}
    cluster_by_seq = collections.defaultdict(list)
    for i in range(len(clusters)):
        seq_id = clusters.sequence_id[i]
        cluster_id = clusters.cluster_id[i]
        cluster_by_seq[seq_id].append(
            (
                clusters.start[i],
                clusters.end[i],
                clusters.cluster_id[i],
            )
        )
        if not "type" in clusters.data.columns:
            cluster_types[cluster_id] = None
        elif clusters.type[i] == "Unknown" or clusters.type[i] is None:
            cluster_types[cluster_id] = ClusterType()
        else:
            cluster_types[cluster_id] = ClusterType(*clusters.type[i].split(";"))

    logger.info("Extracting", "genes belonging to clusters")
    genes_by_cluster = collections.defaultdict(list)
    for seq_id, seq_genes in itertools.groupby(
        genes, key=operator.attrgetter("source.id")
    ):
        for gene in seq_genes:
            for cluster_start, cluster_end, cluster_id in cluster_by_seq[seq_id]:
                if cluster_start <= gene.end and gene.start <= cluster_end:
                    genes_by_cluster[cluster_id].append(gene)

    return [
        Cluster(cluster_id, genes_by_cluster[cluster_id], cluster_types[cluster_id])
        for cluster_id in sorted(clusters.cluster_id)
        if genes_by_cluster[cluster_id]
    ]


def _save_domain_compositions(
    logger: ConsoleLogger,
    all_possible: List[str],
    clusters: List["Cluster"],
    *,
    output_dir: pathlib.Path,
) -> None:
    import numpy
    import scipy.sparse

    logger.info("Saving", "training matrix labels for type classifier")
    with open(os.path.join(output_dir, "domains.tsv"), "w") as out:
        out.writelines(f"{domain}\n" for domain in all_possible)
    with open(os.path.join(output_dir, "types.tsv"), "w") as out:
        writer = csv.writer(out, dialect="excel-tab")
        for cluster in clusters:
            types = ";".join(sorted(cluster.type.names))
            writer.writerow([cluster.id, types])

    logger.info("Building", "new domain composition matrix")
    comp = numpy.array([c.domain_composition(all_possible) for c in clusters])

    comp_out = os.path.join(output_dir, "compositions.npz")
    logger.info("Saving", "new domain composition matrix to file", repr(comp_out))
    scipy.sparse.save_npz(comp_out, scipy.sparse.coo_matrix(comp))


def run(
    args: argparse.Namespace,
    console: Console,
    crf_type: Type["ClusterCRF"],
    classifier_type: Type["TypeClassifier"],
    default_hmms: Callable[[], Iterable["HMM"]],
) -> int:
    logger = ConsoleLogger(console, quiet=args.quiet, verbose=args.verbose)

    # seed RNG
    _common.seed_rng(logger, args.seed)

    # attempt to create the output directory
    outputs = [
        "model.pkl",
        "model.pkl.md5",
        "domains.tsv",
        "types.tsv",
        "compositions.npz",
    ]
    _common.make_output_directory(logger, args.output_dir, outputs)

    # load features
    genes = list(_common.load_genes(logger, args.genes))
    features = _common.load_features(logger, args.features)

    # label genes
    genes = _common.annotate_genes(logger, genes, features)

    # Sort genes
    logger.info("Sorting", "genes by coordinates", level=2)
    genes.sort(key=operator.attrgetter("source.id", "start", "end"))
    for gene in genes:
        gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

    # filter domains by p-value and/or e-value
    genes = _common.filter_domains(
        logger, genes, e_filter=args.e_filter, p_filter=args.p_filter
    )

    # load clusters and label genes inside clusters
    clusters = _common.load_clusters(logger, args.clusters)
    genes = _common.label_genes(logger, genes, clusters)

    # fit CRF
    crf = _common.fit_model(
        logger,
        genes,
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

    # save model
    logger.info("Saving", f"CRF model to {str(args.output_dir)!r}")
    crf.save(args.output_dir)
    _save_transitions(logger, crf, output_dir=args.output_dir)
    _save_weights(logger, crf, output_dir=args.output_dir)

    # extract the domains names
    logger.info("Finding", "the array of possible protein domains", level=2)
    if crf.significant_features is not None:
        all_possible = sorted(crf.significant_features)
    else:
        all_possible = sorted({d.name for g in genes for d in g.protein.domains})

    # compute domain compositions
    _save_domain_compositions(
        logger,
        all_possible,
        _assign_clusters(logger, genes, clusters),
        output_dir=args.output_dir,
    )
    logger.success("Finished", "training new CRF model", level=0)

    # exit succesfully
    return 0
