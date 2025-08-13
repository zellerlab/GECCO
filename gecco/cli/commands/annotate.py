"""Implementation of the ``gecco annotate`` subcommand."""

import argparse
import os
import typing
from typing import (
    Any,
    BinaryIO,
    Container,
    Collection,
    Dict,
    Iterable,
    Iterator,
    Union,
    Optional,
    List,
    TextIO,
    Mapping,
    Callable,
    Iterable,
    Type,
)

from rich.console import Console

from ..._meta import zopen
from .._log import ConsoleLogger
from . import _parser, _common

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord
    from ...hmmer import HMM
    from ...model import Gene
    from ...orf import ORFFinder


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
    _parser.configure_group_table_output(parser, defaults=defaults)

    parser.set_defaults(run=run)


def run(
    args: argparse.Namespace,
    logger: ConsoleLogger,
    crf_type: Type["ClusterCRF"],
    classifier_type: Type["TypeClassifier"],
    default_hmms: Callable[[], Iterable["HMM"]],
) -> int:  # noqa: D102
    # attempt to create the output directory, checking it doesn't
    # already contain output files (or raise a warning)
    base, _ = os.path.splitext(os.path.basename(args.genome))
    outputs = [f"{base}.features.tsv", f"{base}.genes.tsv"]
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
            _common.write_feature_table([])
        logger.warn("No genes were found")
        return 0

    # annotate genes with protein domains and write results
    genes = _common.annotate_domains(
        logger,
        genes,
        hmm_paths=args.hmms,
        default_hmms=default_hmms(),
        disentangle=args.disentangle,
        jobs=args.jobs,
        bit_cutoffs=args.bit_cutoffs,
        e_filter=args.e_filter,
        p_filter=args.p_filter,
    )
    _common.write_feature_table(
        logger,
        genes,
        genome=args.genome,
        output_dir=args.output_dir,
    )

    # report number of domains found
    ndoms = sum(len(gene.protein.domains) for gene in genes)
    if ndoms:
        logger.success("Found", ndoms, "protein domains", level=0)
    else:
        logger.warn("No protein domains were found")

    # exit successfully
    return 0
