"""Implementation of the ``gecco annotate`` subcommand.
"""

import contextlib
import errno
import glob
import gzip
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
from typing import Any, Dict, Union, Optional, List, TextIO, Mapping

from ._base import Command, CommandExit, InvalidArgument
from .._utils import (
    guess_sequences_format,
    in_context,
    patch_showwarnings,
    ProgressReader,
)


class Annotate(Command):  # noqa: D101

    summary = "annotate protein features of one or several contigs."

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
        return f"""
        gecco annotate - {cls.summary}

        Usage:
            gecco annotate --genome <file> [--hmm <hmm>]... [options]

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


        Parameters - Gene Calling:
            -M, --mask                    Enable unknown region masking to
                                          prevent genes from stretching across
                                          unknown nucleotides.

        Parameters - Domain Annotation:
            -e <e>, --e-filter <e>        the e-value cutoff for protein domains
                                          to be included. This is not stable
                                          across versions, so consider using
                                          a p-value filter instead.
            -p <p>, --p-filter <p>        the p-value cutoff for protein domains
                                          to be included. [default: 1e-9]

        Parameters - Debug:
            --hmm <hmm>                   the path to one or more alternative
                                          HMM file to use (in HMMER format).
        """

    def _check(self) -> typing.Optional[int]:
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
            self.format = self._check_flag("--format", optional=True)
            self.genome = self._check_flag("--genome")
            self.hmm = self._check_flag("--hmm", optional=True)
            self.output_dir = self._check_flag("--output-dir")
            self.mask = self._check_flag("--mask", bool)
        except InvalidArgument:
            raise CommandExit(1)

    def _custom_hmms(self):
        from ...hmmer import HMM

        for path in self.hmm:
            base = os.path.basename(path)
            file = open(path, "rb")
            if base.endswith(".gz"):
                base, _ = os.path.splitext(base)
                file = gzip.GzipFile(fileobj=file)
            base, _ = os.path.splitext(base)
            yield HMM(
                id=base,
                version="?",
                url="?",
                path=path,
                size=1,
                relabel_with=r"s/([^\.]*)(\..*)?/\1/"
            )

    # ---

    _OUTPUT_FILES = ["features.tsv", "genes.tsv"]

    def _make_output_directory(self, extensions: List[str]) -> None:
        # Make output directory
        self.info("Using", "output folder", repr(self.output_dir), level=1)
        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except OSError as err:
            self.error("Could not create output directory: {}", err)
            raise CommandExit(err.errno) from err

        # Check if output files already exist
        base, _ = os.path.splitext(os.path.basename(self.genome))
        # output_exts = ["features.tsv", "genes.tsv", "clusters.tsv"]
        # if self.antismash_sideload:
        #     output_exts.append("sideload.json")
        for ext in extensions:
            if os.path.isfile(os.path.join(self.output_dir, f"{base}.{ext}")):
                self.warn("Output folder contains files that will be overwritten")
                break

    def _load_sequences(self):
        from Bio import SeqIO

        try:
            # guess format or use the one given in CLI
            if self.format is not None:
                format = self.format
                self.info("Using", "user-provided sequence format", repr(format), level=2)
            else:
                self.info("Detecting", "sequence format from file contents", level=2)
                format = guess_sequences_format(self.genome)
                self.success("Detected", "format of input as", repr(format), level=2)
            # get filesize and unit
            input_size = os.stat(self.genome).st_size
            total, scale, unit = ProgressReader.scale_size(input_size)
            task = self.progress.add_task("Loading sequences", total=total, unit=unit, precision=".1f")
            # load sequences
            self.info("Loading", "sequences from genomic file", repr(self.genome), level=1)
            with ProgressReader(open(self.genome, "rb"), self.progress, task, scale) as f:
                sequences = list(SeqIO.parse(io.TextIOWrapper(f), format))
        except FileNotFoundError as err:
            self.error("Could not find input file:", repr(self.genome))
            raise CommandExit(err.errno) from err
        except ValueError as err:
            self.error("Failed to load sequences:", err)
            raise CommandExit(getattr(err, "errno", 1)) from err
        else:
            self.success("Found", len(sequences), "sequences", level=1)
            return sequences

    def _extract_genes(self, sequences):
        from ...orf import PyrodigalFinder

        self.info("Extracting", "genes from input sequences", level=1)
        orf_finder = PyrodigalFinder(metagenome=True, mask=self.mask, cpus=self.jobs)

        unit = "contigs" if len(sequences) > 1 else "contig"
        task = self.progress.add_task(description="Finding ORFs", total=len(sequences), unit=unit, precision="")

        def callback(record, found, total):
            self.success("Found", found, "genes in record", repr(record.id), level=2)
            self.progress.update(task, advance=1)

        return list(orf_finder.find_genes(sequences, progress=callback))

    def _annotate_domains(self, genes, whitelist=None):
        from ...hmmer import PyHMMER, embedded_hmms

        self.info("Running", "HMMER domain annotation", level=1)

        # Run all HMMs over ORFs to annotate with protein domains
        hmms = list(self._custom_hmms() if self.hmm else embedded_hmms())
        task = self.progress.add_task(description=f"Annotating domains", unit="HMMs", total=len(hmms), precision="")
        for hmm in self.progress.track(hmms, task_id=task, total=len(hmms)):
            task = self.progress.add_task(description=f"  {hmm.id} v{hmm.version}", total=hmm.size, unit="domains", precision="")
            callback = lambda h, t: self.progress.update(task, advance=1)
            self.info("Starting", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
            PyHMMER(hmm, self.jobs, whitelist).run(genes, progress=callback)
            self.success("Finished", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
            self.progress.update(task_id=task, visible=False)

        # Count number of annotated domains
        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.success("Found", count, "domains across all proteins", level=1)

        # Filter i-evalue and p-value if required
        genes = self._filter_domains(genes)

        # Sort genes
        self.info("Sorting", "genes by coordinates", level=2)
        genes.sort(key=lambda g: (g.source.id, g.start, g.end))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

        return genes

    def _filter_domains(self, genes):
        # Filter i-evalue and p-value if required
        if self.e_filter is not None:
            self.info("Excluding", "domains with e-value over", self.e_filter, level=1)
            key = lambda d: d.i_evalue < self.e_filter
            genes = [
                gene.with_protein(gene.protein.with_domains(filter(key, gene.protein.domains)))
                for gene in genes
            ]
        if self.p_filter is not None:
            self.info("Excluding", "domains with p-value over", self.p_filter, level=1)
            key = lambda d: d.pvalue < self.p_filter
            genes = [
                gene.with_protein(gene.protein.with_domains(filter(key, gene.protein.domains)))
                for gene in genes
            ]
        if self.p_filter is not None or self.e_filter is not None:
            count = sum(1 for gene in genes for domain in gene.protein.domains)
            self.info("Using", "remaining", count, "domains", level=1)
        return genes

    def _write_feature_table(self, genes):
        from ...model import FeatureTable

        base, _ = os.path.splitext(os.path.basename(self.genome))
        pred_out = os.path.join(self.output_dir, f"{base}.features.tsv")
        self.info("Writing", "feature table to", repr(pred_out), level=1)
        with open(pred_out, "w") as f:
            FeatureTable.from_genes(genes).dump(f)

    def _write_genes_table(self, genes):
        from ...model import GeneTable

        base, _ = os.path.splitext(os.path.basename(self.genome))
        pred_out = os.path.join(self.output_dir, f"{base}.genes.tsv")
        self.info("Writing", "gene table to", repr(pred_out), level=1)
        with open(pred_out, "w") as f:
            GeneTable.from_genes(genes).dump(f)

    # ---

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check the CLI arguments were fine and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # attempt to create the output directory, checking it doesn't
            # already contain output files (or raise a warning)
            self._make_output_directory(extensions=["features.tsv", "genes.tsv"])
            # load sequences and extract genes
            sequences = self._load_sequences()
            genes = self._extract_genes(sequences)
            self._write_genes_table(genes)
            if genes:
                self.success("Found", "a total of", len(genes), "genes", level=1)
            else:
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
