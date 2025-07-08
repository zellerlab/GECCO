"""Implementation of the ``gecco run`` subcommand."""

import argparse
import os
import typing
import pathlib
from typing import Optional

from rich.console import Console

from ..._meta import zopen
from .._log import ConsoleLogger
from . import _parser, _common


def configure_parser(parser: argparse.ArgumentParser):
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit.",
    )

    params_arguments = parser.add_argument_group(
        "Arguments",
    )
    params_arguments.add_argument(
        "-g",
        "--genome",
        required=True,
        help=(
            "A genomic file containing one or more sequences to use as input. "
            "Must be in one of the sequence formats supported by Biopython"
        ),
    )

    params_parameters = parser.add_argument_group(
        "Parameters",
    )
    params_parameters.add_argument(
        "-f",
        "--format",
        help=(
            "The format of the input file, as a Biopython format string. "
            "GECCO is able to recognize FASTA and GenBank files automatically "
            "if this is not given."
        ),
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

    params_gene_calling = _parser.configure_group_gene_calling(parser)
    params_domain_annotation = _parser.configure_group_domain_annotation(parser)
    params_cluster_detection = _parser.configure_group_cluster_detection(parser)

    params_output = _parser.configure_group_table_output(parser)
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
        help=("Write an AntiSMASH v6 sideload JSON file next to the output " "files."),
    )

    parser.set_defaults(run=run)


def run(args: argparse.Namespace, console: Console) -> int:
    logger = ConsoleLogger(console, quiet=args.quiet, verbose=args.verbose)

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

    # generate whitelist from internal feature list of model
    whitelist = _common.load_model_domains(logger, model=args.model)

    # annotate genes with protein domains
    genes = _common.annotate_domains(
        logger,
        genes,
        hmm_paths=args.hmms,
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

    # predict types for putative clusters
    classifier = _common.load_type_classifier(logger, model=args.model)
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
    logger.success("Found", len(clusters), "biosynthetic gene", unit, level=0)

    # exit successfully
    return 0
