"""Implementation of the ``gecco run`` subcommand."""

import argparse
import os
import typing
import pathlib
from typing import Optional, Type, Iterable, Callable, Dict

from rich.console import Console

from ..._meta import zopen
from .._log import ConsoleLogger
from . import _parser, _common


def configure_parser(
    parser: argparse.ArgumentParser, 
    console: Console,
    program: str,
    version: str,
    *,
    defaults: Dict[str, object],
):
    _parser.configure_common(parser, console, program, version, defaults=defaults)

    _parser.configure_group_input_sequences(parser, defaults=defaults)
    _parser.configure_group_gene_calling(parser, defaults=defaults)
    _parser.configure_group_domain_annotation(parser, defaults=defaults)
    _parser.configure_group_cluster_detection(parser, defaults=defaults)

    params_output = _parser.configure_group_table_output(parser, defaults=defaults)
    params_output.add_argument(
        "--merge-gbk",
        action="store_true",
        help=(
            "Output a single GenBank file containing every detected cluster "
            "instead of writing one file per cluster."
        ),
    )
    params_output.add_argument(
        "--antismash-sideload",
        action="store_true",
        help=("Write an AntiSMASH v6 sideload JSON file next to the output files."),
    )

    parser.set_defaults(run=run)


def run(
    args: argparse.Namespace,
    logger: ConsoleLogger,
    crf_type: Type["ClusterCRF"],
    classifier_type: Type["TypeClassifier"],
    default_hmms: Callable[[], Iterable["HMM"]],
) -> int:
    # attempt to create the output directory, checking it doesn't
    # already contain output files (or raise a warning)
    base, _ = os.path.splitext(os.path.basename(args.genome))
    outputs = [f"{base}.features.tsv", f"{base}.genes.tsv", f"{base}.clusters.tsv"]
    if args.antismash_sideload:
        outputs.append(f"{base}.sideload.json")
    if args.merge_gbk:
        outputs.append(f"{base}.clusters.gbk")
    _common.make_output_directory(logger, args.output_dir, outputs)

    # load sequences and extract genes
    sequences = list(
        _common.load_sequences(
            logger,
            args.genome,
            format=args.format,
        )
    )
    genes = list(
        _common.extract_genes(
            logger,
            sequences,
            cds_feature=args.cds_feature,
            locus_tag=args.locus_tag,
            mask=args.mask,
            jobs=args.jobs,
        )
    )

    # write gene table
    _common.write_genes_table(
        logger, genes, genome=args.genome, output_dir=args.output_dir
    )
    if genes:
        logger.success("Found", "a total of", len(genes), "genes", level=1)
    else:
        if args.force_tsv:
            _common.write_gene_table(
                logger, [], genome=args.genome, output_dir=args.output_dir
            )
            _common.write_feature_table(
                logger, [], genome=args.genome, output_dir=args.output_dir
            )
            _common.write_cluster_table(
                logger, [], genome=args.genome, output_dir=args.output_dir
            )
        logger.warn("No genes were found")
        return 0

    # load type classifier (to access compositional data)
    classifier = _common.load_type_classifier(
        logger,
        model=args.model,
        classifier_type=classifier_type,
    )

    # generate whitelist from internal feature list of model
    whitelist = _common.load_model_domains(
        logger,
        classifier,
    )

    # annotate genes with protein domains
    genes = _common.annotate_domains(
        logger,
        genes,
        hmm_paths=args.hmms,
        default_hmms=default_hmms(),
        whitelist=whitelist,
        disentangle=args.disentangle,
        jobs=args.jobs,
        bit_cutoffs=args.bit_cutoffs,
        e_filter=args.e_filter,
        p_filter=args.p_filter,
    )

    # predict gene cluster probabilities
    genes = _common.predict_probabilities(
        logger,
        genes,
        model=args.model,
        pad=args.pad,
        crf_type=crf_type,
    )

    # write gene and feature tables
    _common.write_genes_table(
        logger,
        genes,
        genome=args.genome,
        output_dir=args.output_dir,
    )
    _common.write_feature_table(
        logger,
        genes,
        genome=args.genome,
        output_dir=args.output_dir,
    )

    # extract clusters from probability vector
    clusters = _common.extract_clusters(
        logger,
        genes,
        threshold=args.threshold,
        postproc=args.postproc,
        cds=args.cds,
        edge_distance=args.edge_distance,
    )
    if clusters:
        logger.success("Found", len(clusters), "potential gene clusters", level=1)
    else:
        logger.warn("No gene clusters were found")
        if args.force_tsv:
            _common.write_cluster_table(
                logger,
                clusters,
                genome=args.genome,
                output_dir=args.output_dir,
            )
        return 0

    # predict types for putative clusters if the type classifier has more
    # than one class to predict from
    if len(classifier.classes_) > 1:
        clusters = _common.predict_types(logger, clusters, classifier=classifier)

    # write results
    logger.info(
        "Writing", "result files to folder", repr(str(args.output_dir)), level=1
    )
    _common.write_cluster_table(
        logger, clusters, genome=args.genome, output_dir=args.output_dir
    )
    _common.write_clusters(
        logger,
        clusters,
        merge=args.merge_gbk,
        genome=args.genome,
        output_dir=args.output_dir,
    )
    if args.antismash_sideload:
        _common.write_sideload_json(clusters)
    unit = "cluster" if len(clusters) == 1 else "clusters"
    logger.success("Found", len(clusters), "gene", unit, level=0)

    # exit successfully
    return 0
