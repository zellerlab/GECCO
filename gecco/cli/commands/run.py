"""Implementation of the ``gecco run`` subcommand.
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

import numpy
import rich.emoji
import rich.progress
from Bio import SeqIO

from ._base import Command
from .._utils import guess_sequences_format, in_context, patch_showwarnings
from ...crf import ClusterCRF
from ...hmmer import PyHMMER, HMMER, HMM, embedded_hmms
from ...model import FeatureTable, ClusterTable, ProductType
from ...orf import PyrodigalFinder
from ...types import TypeClassifier
from ...refine import ClusterRefiner


class Run(Command):  # noqa: D101

    summary = "predict BGC from one or several contigs."
    doc = f"""
    gecco run - {summary}

    Usage:
        gecco run (-h | --help)
        gecco run --genome <file> [--hmm <hmm>]... [options]

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
        -o <out>, --output-dir <out>  the directory in which to write the
                                      output files. [default: .]
        -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                      multithreading. Use 0 to use all of the
                                      available CPUs. [default: 0]

    Parameters - Domain Annotation:
        -e <e>, --e-filter <e>        the e-value cutoff for protein domains
                                      to be included. [default: 1e-5]

    Parameters - Cluster Detection:
        -c, --cds <N>                 the minimum number of coding sequences a
                                      valid cluster must contain. [default: 3]
        -m <m>, --threshold <m>       the probability threshold for cluster
                                      detection. Default depends on the
                                      post-processing method (0.4 for gecco,
                                      0.6 for antismash).
        --postproc <method>           the method to use for cluster validation
                                      (antismash or gecco). [default: gecco]

    Parameters - Debug:
        --model <directory>           the path to an alternative CRF model
                                      to use (obtained with `gecco train`).
        --hmm <hmm>                   the path to one or more alternative HMM
                                      file to use (in HMMER format).
    """

    def __init__(
        self,
        argv: Optional[List[str]] = None,
        stream: Optional[TextIO] = None,
        logger: Optional[logging.Logger] = None,
        options: Optional[Mapping[str, Any]] = None,
        config: Optional[Dict[Any, Any]] = None,
    ) -> None:
        super().__init__(argv, stream, logger, options, config)
        self.progress = rich.progress.Progress(
            rich.progress.SpinnerColumn(finished_text="[green]:heavy_check_mark:[/]"),
            "[progress.description]{task.description}",
            rich.progress.BarColumn(bar_width=60),
            "[progress.completed]{task.completed}/{task.total}",
            "[progress.completed]{task.fields[unit]}",
            "[progress.percentage]{task.percentage:>3.0f}%",
            rich.progress.TimeElapsedColumn(),
            rich.progress.TimeRemainingColumn(),
            console=self.console,
            disable=self.quiet > 0,
        )
        self.console = self.progress.console

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check value of numeric arguments
        self.args["--cds"] = int(self.args["--cds"])
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for --e-filter`: {}", e_filter)
            return 1

        # Use default threshold value dependeing on postprocessing method
        if self.args["--threshold"] is None:
            if self.args["--postproc"] == "gecco":
                self.args["--threshold"] = 0.4
            elif self.args["--postproc"] == "antismash":
                self.args["--threshold"] = 0.6
        else:
            self.args["--threshold"] = float(self.args["--threshold"])

        # Check the `--cpu`flag
        try:
            self.jobs = int(self.args["--jobs"])
        except ValueError:
            self.error("Invalid value for --jobs:", repr(self.args["--jobs"]))

        return None

    def _custom_hmms(self):
        for path in self.args["--hmm"]:
            base = os.path.basename(path)
            if base.endswith(".gz"):
                base, _ = os.path.splitext(base)
            base, _ = os.path.splitext(base)
            yield HMM(
                id=base,
                version="?",
                url="?",
                path=path,
                relabel_with=r"s/([^\.]*)(\..*)?/\1/"
            )

    # ---

    def make_output_directory(self):
        # Make output directory
        out_dir = self.args["--output-dir"]

        self.info("Using", "output folder", repr(out_dir), level=1)
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as e:
            self.error("Could not create output directory: {}", e)
            return e.errno

        # Check if output files already exist
        base, _ = os.path.splitext(os.path.basename(self.args["--genome"]))
        for ext in ["features.tsv", "clusters.tsv"]:
            if os.path.isfile(os.path.join(out_dir, f"{base}.{ext}")):
                self.warn("Output folder contains files that will be overwritten")
                break

    def load_sequences(self):
        genome = self.args["--genome"]

        if self.args["--format"] is not None:
            format = self.args["--format"]
            self.info("Using", "user-provided sequence format", repr(format), level=2)
        else:
            self.info("Detecting", "sequence format from file contents", level=2)
            format = guess_sequences_format(genome)
            self.success("Detected", "format of input as", repr(format), level=2)

        self.info("Loading", "sequences from genomic file", repr(genome), level=1)
        try:
            sequences = list(SeqIO.parse(genome, format))
        except FileNotFoundError:
            self.error()
        except Exception as err:
            self.error("Failed to load sequences:", err)
            return getattr(err, "errno", 1)
        else:
            self.success("Found", len(sequences), "sequences", level=1)
            return sequences

    def extract_genes(self, sequences):
        self.info("Extracting", "genes from input sequences", level=1)
        orf_finder = PyrodigalFinder(metagenome=True, cpus=int(self.args["--jobs"]))

        unit = "contigs" if len(sequences) > 1 else "contig"
        task = self.progress.add_task(description="ORFs finding", total=len(sequences), unit=unit)

        def callback(record, found, total):
            self.success("Found", found, "genes in record", repr(record.id), level=2)
            self.progress.update(task, advance=1)

        return list(orf_finder.find_genes(sequences, progress=callback))

    def annotate_domains(self, genes):
        self.info("Running", "HMMER domain annotation", level=1)

        # Run all HMMs over ORFs to annotate with protein domains
        hmms = list(self._custom_hmms() if self.args["--hmm"] else embedded_hmms())
        task = self.progress.add_task(description=f"HMM annotation", unit="HMMs", total=len(hmms))
        for hmm in self.progress.track(hmms, task_id=task):
            task = self.progress.add_task(description=f"{hmm.id} v{hmm.version}", total=1, unit="domains")
            callback = lambda n, total: self.progress.update(task, advance=1, total=total)
            self.info("Starting", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
            features = PyHMMER(hmm, int(self.args["--jobs"])).run(genes, progress=callback)
            self.success("Finished", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
            self.progress.update(task_id=task, visible=False)

        # Count number of annotated domains
        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.success("Found", count, "domains across all proteins", level=1)

        # Filter i-evalue
        self.info("Filtering", "results with e-value under", self.args["--e-filter"], level=1)
        for gene in genes:
            key = lambda d: d.i_evalue < self.args["--e-filter"]
            gene.protein.domains = list(filter(key, gene.protein.domains))

        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.info("Using", "remaining", count, "domains", level=2)

        # Sort genes
        self.info("Sorting", "genes by coordinates", level=2)
        genes.sort(key=lambda g: (g.source.id, g.start, g.end))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

        return genes

    def predict_probabilities(self, genes):
        if self.args["--model"] is None:
            self.info("Loading", "embedded CRF pre-trained model", level=1)
        else:
            self.info("Loading", "CRF pre-trained model from", repr(self.args["--model"]), level=1)
        crf = ClusterCRF.trained(self.args["--model"])

        self.info("Predicting", "cluster probabilitites with the CRF model", level=1)
        unit = "genes" if len(genes) > 1 else "gene"
        task = self.progress.add_task("Prediction", total=len(genes), unit=unit)
        return list(crf.predict_probabilities(self.progress.track(genes, task_id=task)))

    def write_feature_table(self, genes):
        base, _ = os.path.splitext(os.path.basename(self.args["--genome"]))
        pred_out = os.path.join(self.args["--output-dir"], f"{base}.features.tsv")
        self.info("Writing", "feature table to", repr(pred_out), level=1)
        with open(pred_out, "w") as f:
            FeatureTable.from_genes(genes).dump(f)

    def extract_clusters(self, genes):
        self.info("Extracting", "predicted biosynthetic regions", level=1)
        refiner = ClusterRefiner(self.args["--threshold"], self.args["--postproc"], self.args["--cds"])

        total = len({gene.source.id for gene in genes})
        unit = "contigs" if total > 1 else "contig"
        task = self.progress.add_task("Segmentation", total=total, unit=unit)

        clusters = []
        gene_groups = itertools.groupby(genes, lambda g: g.source.id)
        for _, gene_group in self.progress.track(gene_groups, task_id=task, total=total):
            clusters.extend(refiner.iter_clusters(gene_group))

        return clusters

    def predict_types(self, clusters):
        self.info("Predicting", "BGC types", level=1)

        unit = "cluster" if len(clusters) == 1 else "clusters"
        task = self.progress.add_task("Type prediction", total=len(clusters), unit=unit)

        clusters_new = []
        classifier = TypeClassifier.trained(self.args["--model"])
        for cluster in self.progress.track(clusters, task_id=task):
            clusters_new.extend(classifier.predict_types([cluster]))
            if cluster.type != ProductType.Unknown:
                name = "/".join(ty.name for ty in cluster.type.unpack())
                prob = "/".join(f"[purple]{cluster.type_probabilities[ty]:.0%}[/]" for ty in cluster.type.unpack())
                self.success(f"Predicted type of [bold blue]{cluster.id}[/] as [bold blue]{name}[/] ({prob} confidence)")
            else:
                ty = max(cluster.type_probabilities, key=cluster.type_probabilities.get)
                prob = f"[purple]{cluster.type_probabilities[ty]}[/]"
                self.warn(f"Couldn't assign type to [bold blue]{cluster.id}[/] (maybe [bold blue]{ty}[/], {prob} confidence")

        return clusters_new

    def write_cluster_table(self, clusters):
        base, _ = os.path.splitext(os.path.basename(self.args["--genome"]))
        cluster_out = os.path.join(self.args["--output-dir"], f"{base}.clusters.tsv")
        self.info("Writing", "cluster table to", repr(cluster_out), level=1)
        with open(cluster_out, "w") as out:
            ClusterTable.from_clusters(clusters).dump(out)

    def write_clusters(self, clusters):
        for cluster in clusters:
            gbk_out = os.path.join(self.args["--output-dir"], f"{cluster.id}.gbk")
            self.info("Writing", f"cluster [bold blue]{cluster.id}[/] to", repr(gbk_out), level=1)
            SeqIO.write(cluster.to_seq_record(), gbk_out, "genbank")

    # ---

    @in_context
    def __call__(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        retcode = self._check()
        if retcode is not None:
            return retcode

        ctx.enter_context(self.progress)
        ctx.enter_context(patch_showwarnings(self._showwarnings))

        retcode = self.make_output_directory()
        if retcode:
            return retcode

        sequences = self.load_sequences()
        if not isinstance(sequences, list):
            return sequences

        genes = self.extract_genes(sequences)
        if genes:
            self.success("Found", "a total of", len(genes), "genes", level=1)
        else:
            self.warn("No genes were found")
            return 0

        genes = self.annotate_domains(genes)

        genes = self.predict_probabilities(genes)
        self.write_feature_table(genes)

        clusters = self.extract_clusters(genes)
        if clusters:
            self.success("Found", len(clusters), "potential gene clusters", level=1)
        else:
            self.warn("No gene clusters were found")
            return 0

        clusters = self.predict_types(clusters)

        self.info("Writing", "result files to folder", repr(self.args["--output-dir"]), level=1)
        self.write_cluster_table(clusters)
        self.write_clusters(clusters)

        self.success("Found", len(clusters), "biosynthetic gene clusters", level=0)
        return 0
