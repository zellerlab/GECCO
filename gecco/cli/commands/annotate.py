"""Implementation of the ``gecco annotate`` subcommand.
"""

import contextlib
import errno
import glob
import itertools
import io
import logging
import multiprocessing
import operator
import os
import pickle
import tempfile
import typing
import signal
from typing import Any, BinaryIO, Container, Dict, Iterable, Union, Optional, List, TextIO, Mapping

from ._base import Command, CommandExit, InvalidArgument
from ._mixins import SequenceLoaderMixin, OutputWriterMixin, AnnotatorMixin
from .._utils import (
    guess_sequences_format,
    in_context,
    patch_showwarnings,
    ProgressReader,
)

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord
    from ...hmmer import HMM
    from ...model import Gene
    from ...orf import ORFFinder


class Annotate(SequenceLoaderMixin, OutputWriterMixin, AnnotatorMixin):  # noqa: D101

    summary = "annotate protein features of one or several contigs."

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        return f"""
        gecco annotate - {cls.summary}

        Usage:
            gecco annotate --genome <file> [--hmm <hmm>]... [--hmm-x <hmm>]... [options]

        Arguments:
            -g <file>, --genome <file>    a genomic file containing one or more
                                          sequences to use as input. Must be in
                                          one of the sequences format supported
                                          by Biopython.

        Parameters:
            -f <fmt>, --format <fmt>      the format of the input file, as a
                                          Biopython format string. GECCO is able
                                          to recognize FASTA and GenBank files
                                          automatically if this is not given.
            -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                          multithreading. Use 0 to use all of
                                          the available CPUs. [default: 0]

        Parameters - Output:
            -o <out>, --output-dir <out>  the directory in which to write the
                                          output files. [default: .]
            --force-tsv                   always write TSV output files even
                                          when they are empty (e.g. because
                                          no genes or no clusters were found).


        Parameters - Gene Calling:
            -M, --mask                    Enable unknown region masking to
                                          prevent genes from stretching across
                                          unknown nucleotides.
            --cds-feature <cds_feature>   Extract genes from annotated records
                                          using a feature rather than calling
                                          genes from scratch.
            --locus-tag <locus_tag>       The name of the feature qualifier
                                          to use for naming extracted genes
                                          when using the ``--cds-feature``
                                          flag. [default: locus_tag]

        Parameters - Domain Annotation:
            --hmm <hmm>                   the path to one or more alternative
                                          HMM file to use (in HMMER format).
            -e <e>, --e-filter <e>        the e-value cutoff for protein domains
                                          to be included. This is not stable
                                          across versions, so consider using
                                          a p-value filter instead.
            -p <p>, --p-filter <p>        the p-value cutoff for protein domains
                                          to be included. [default: 1e-9]
            --bit-cutoffs <name>          use bitscore cutoffs (one of *noise*,
                                          *gathering*, or *trusted*) to filter
                                          domain annotations.
            --disentangle                 disentangle overlapping domains in 
                                          each gene by keeping only the domains
                                          with the lowest E-value over a given
                                          position.

        """

    def _check(self) -> None:
        _BITSCORE_CUTOFFS = {"gathering", "noise", "trusted"}
        super()._check()
        try:
            self.e_filter = self._check_flag(
                "--e-filter",
                float,
                lambda x: x > 0,
                hint="real number above 0",
                optional=True,
            )
            self.p_filter = self._check_flag(
                "--p-filter",
                float,
                lambda x: x > 0,
                hint="real number above 0",
                optional=True,
            )
            self.jobs = self._check_flag("--jobs", int, lambda x: x >= 0, hint="positive or null integer")
            self.format: Optional[str] = self._check_flag("--format", optional=True)
            self.genome: str = self._check_flag("--genome")
            self.hmm: Optional[List[str]] = self._check_flag("--hmm", optional=True)
            self.output_dir: str = self._check_flag("--output-dir")
            self.mask = self._check_flag("--mask", bool)
            self.force_tsv = self._check_flag("--force-tsv", bool)
            self.cds_feature: Optional[str] = self._check_flag("--cds-feature", optional=True)
            self.locus_tag: str = self._check_flag("--locus-tag")
            self.disentangle = self._check_flag("--disentangle", bool)
            self.bit_cutoffs: str = self._check_flag(
                "--bit-cutoffs",
                str,
                _BITSCORE_CUTOFFS.__contains__,
                optional=True,
                hint="one of {}".format(", ".join(sorted(_BITSCORE_CUTOFFS)))
            )
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    _OUTPUT_FILES = ["features.tsv", "genes.tsv"]

    def _extract_genes(self, sequences: List["SeqRecord"]) -> List["Gene"]:
        from ...orf import PyrodigalFinder, CDSFinder

        self.info("Extracting", "genes from input sequences", level=1)
        if self.cds_feature is None:
            self.info("Using", "Pyrodigal in metagenomic mode", level=2)
            orf_finder: ORFFinder = PyrodigalFinder(metagenome=True, mask=self.mask, cpus=self.jobs)
        else:
            self.info("Using", f"record features named {self.cds_feature!r}", level=2)
            orf_finder = CDSFinder(feature=self.cds_feature, locus_tag=self.locus_tag)

        unit = "contigs" if len(sequences) > 1 else "contig"
        task = self.progress.add_task(description="Finding ORFs", total=len(sequences), unit=unit, precision="")

        def callback(record: "SeqRecord", found: int) -> None:
            self.success("Found", found, "genes in record", repr(record.id), level=2)
            self.progress.update(task, advance=1)

        return list(orf_finder.find_genes(sequences, progress=callback))

    # ---

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check the CLI arguments were fine and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))  # type: ignore
            # attempt to create the output directory, checking it doesn't
            # already contain output files (or raise a warning)
            base, _ = os.path.splitext(os.path.basename(self.genome))
            outputs = [f"{base}.features.tsv", f"{base}.genes.tsv"]
            self._make_output_directory(outputs)
            # load sequences and extract genes
            sequences = list(self._load_sequences())
            genes = self._extract_genes(sequences)
            self._write_genes_table(genes)
            if genes:
                self.success("Found", "a total of", len(genes), "genes", level=1)
            else:
                if self.force_tsv:
                    self._write_feature_table([])
                self.warn("No genes were found")
                return 0
            # annotate domains and write results
            genes = self._annotate_domains(genes)
            self._write_feature_table(genes)
            ndoms = sum(1 for gene in genes for domain in gene.protein.domains)
            # report number of proteins found
            if ndoms:
                self.success("Found", ndoms, "protein domains", level=0)
            else:
                self.warn("No protein domains were found")
        except CommandExit as cexit:
            return cexit.code
        except KeyboardInterrupt:
            self.error("Interrupted")
            return -signal.SIGINT
        except Exception as err:
            self.progress.stop()
            raise
        else:
            return 0
