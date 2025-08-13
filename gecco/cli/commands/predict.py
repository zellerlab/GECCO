import argparse
import os
import operator
from typing import Type, Callable, Iterable, Dict

from rich.console import Console

from . import _parser, _common
from .._log import ConsoleLogger


def configure_parser(
    parser: argparse.ArgumentParser, 
    console: Console,
    program: str,
    version: str,
    *,
    defaults: Dict[str, object],
):
    _parser.configure_common(parser, console, program, version, defaults=defaults)

    _parser.configure_group_input_sequences(parser, short=False, defaults=defaults)
    _parser.configure_group_input_tables(parser, clusters=False, defaults=defaults)
    _parser.configure_group_domain_filter(parser, defaults=defaults)
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
        help=("Write an AntiSMASH v6 sideload JSON file next to the output " "files."),
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
    outputs = [f"{base}.clusters.tsv", f"{base}.features.tsv", f"{base}.genes.tsv"]
    if args.antismash_sideload:
        outputs.append(f"{base}.sideload.json")
    if args.merge_gbk:
        outputs.append(f"{base}.clusters.gbk")
    _common.make_output_directory(logger, args.output_dir, outputs)

    # load features
    genes = list(_common.load_genes(logger, args.genes))
    features = _common.load_features(logger, args.features)

    # load sequences
    sequences = list(
        _common.load_sequences(
            logger,
            args.genome,
            format=args.format,
        )
    )

    # label genes
    genes = _common.annotate_genes(logger, genes, features)
    genes = list(_common.assign_sources(logger, sequences, genes, genome=args.genome))

    # Sort genes
    logger.info("Sorting", "genes by coordinates", level=2)
    genes.sort(key=operator.attrgetter("source.id", "start", "end"))
    for gene in genes:
        gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

    # filter domains by p-value and/or e-value
    genes = _common.filter_domains(
        logger, genes, e_filter=args.e_filter, p_filter=args.p_filter
    )

    # predict probabilities with CRF model
    genes = _common.predict_probabilities(
        logger, genes, model=args.model, pad=args.pad, crf_type=crf_type
    )
    _common.write_genes_table(
        logger, genes, genome=args.genome, output_dir=args.output_dir
    )
    _common.write_feature_table(
        logger, genes, genome=args.genome, output_dir=args.output_dir
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
    classifier = _common.load_type_classifier(
        logger,
        model=args.model,
        classifier_type=classifier_type,
    )
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
