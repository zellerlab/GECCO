"""Implementation of the ``gecco annotate`` subcommand.
"""

import contextlib
import errno
import glob
import itertools
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
from .._utils import guess_sequences_format, in_context, patch_showwarnings


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
            -o <out>, --output <out>      the file where to write the feature
                                          table. [default: features.tsv]
            -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                          multithreading. Use 0 to use all of the
                                          available CPUs. [default: 0]

        Parameters - Domain Annotation:
            -e <e>, --e-filter <e>        the e-value cutoff for protein domains
                                          to be included. [default: 1e-5]

        Parameters - Debug:
            --hmm <hmm>                   the path to one or more alternative
                                          HMM file to use (in HMMER format).
        """

    def _check(self) -> typing.Optional[int]:
        super()._check()
        try:
            self.e_filter = self._check_flag("--e-filter", float, lambda x: 0 <= x <= 1, hint="real number between 0 and 1")
            self.jobs = self._check_flag("--jobs", int, lambda x: x >= 0, hint="positive or null integer")
            self.format = self._check_flag("--format")
            self.genome = self._check_flag("--genome")
            self.hmm = self._check_flag("--hmm")
            self.output = self._check_flag("--output")
        except InvalidArgument:
            raise CommandExit(1)

    def _custom_hmms(self):
        import pyhmmer
        from ...hmmer import HMM

        for path in self.hmm:
            base = os.path.basename(path)
            if base.endswith(".gz"):
                base, _ = os.path.splitext(base)
            base, _ = os.path.splitext(base)
            with pyhmmer.plan7.HMMFile(path) as hmm_file:
                size = sum(1 for _ in hmm_file)
            yield HMM(
                id=base,
                version="?",
                url="?",
                path=path,
                size=size,
                relabel_with=r"s/([^\.]*)(\..*)?/\1/"
            )

    # ---

    def _load_sequences(self):
        from Bio import SeqIO

        if self.format is not None:
            format = self.format
            self.info("Using", "user-provided sequence format", repr(format), level=2)
        else:
            self.info("Detecting", "sequence format from file contents", level=2)
            format = guess_sequences_format(self.genome)
            self.success("Detected", "format of input as", repr(format), level=2)

        self.info("Loading", "sequences from genomic file", repr(self.genome), level=1)
        try:
            sequences = list(SeqIO.parse(self.genome, format))
        except FileNotFoundError as err:
            self.error("Could not find input file:", repr(self.genome))
            raise CommandExit(e.errno) from err
        except Exception as err:
            self.error("Failed to load sequences:", err)
            raise CommandExit(getattr(err, "errno", 1)) from err
        else:
            self.success("Found", len(sequences), "sequences", level=1)
            return sequences

    def _extract_genes(self, sequences):
        from ...orf import PyrodigalFinder

        self.info("Extracting", "genes from input sequences", level=1)
        orf_finder = PyrodigalFinder(metagenome=True, cpus=self.jobs)

        unit = "contigs" if len(sequences) > 1 else "contig"
        task = self.progress.add_task(description="ORFs finding", total=len(sequences), unit=unit)

        def callback(record, found, total):
            self.success("Found", found, "genes in record", repr(record.id), level=2)
            self.progress.update(task, advance=1)

        return list(orf_finder.find_genes(sequences, progress=callback))

    def _annotate_domains(self, genes):
        from ...hmmer import PyHMMER, embedded_hmms

        self.info("Running", "HMMER domain annotation", level=1)

        # Run all HMMs over ORFs to annotate with protein domains
        hmms = list(self._custom_hmms() if self.hmm else embedded_hmms())
        task = self.progress.add_task(description=f"HMM annotation", unit="HMMs", total=len(hmms))
        for hmm in self.progress.track(hmms, task_id=task, total=len(hmms)):
            task = self.progress.add_task(description=f"{hmm.id} v{hmm.version}", total=hmm.size, unit="domains")
            callback = lambda h, t: self.progress.update(task, advance=1)
            self.info("Starting", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
            features = PyHMMER(hmm, self.jobs).run(genes, progress=callback)
            self.success("Finished", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
            self.progress.update(task_id=task, visible=False)

        # Count number of annotated domains
        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.success("Found", count, "domains across all proteins", level=1)

        # Filter i-evalue
        self.info("Filtering", "results with e-value under", self.e_filter, level=1)
        key = lambda d: d.i_evalue < self.e_filter
        for gene in genes:
            gene.protein.domains = list(filter(key, gene.protein.domains))

        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.info("Using", "remaining", count, "domains", level=2)

        # Sort genes
        self.info("Sorting", "genes by coordinates", level=2)
        genes.sort(key=lambda g: (g.source.id, g.start, g.end))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

        return genes

    def _write_feature_table(self, genes):
        from ...model import FeatureTable

        self.info("Writing", "feature table to", repr(self.output), level=1)
        with open(self.output, "w") as f:
            FeatureTable.from_genes(genes).dump(f)

    # ---

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check the CLI arguments were fine and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # load sequences and extract genes
            sequences = self._load_sequences()
            genes = self._extract_genes(sequences)
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
