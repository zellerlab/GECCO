"""Implementation of the ``gecco rerun`` subcommand.
"""

import contextlib
import errno
import glob
import io
import itertools
import json
import logging
import multiprocessing
import operator
import os
import pickle
import tempfile
import typing
import signal
from typing import Any, Dict, Union, Optional, Iterable, List, TextIO, Mapping

import rich.emoji
import rich.progress

from ... import __version__
from ._base import Command, CommandExit, InvalidArgument
from ._mixins import TableLoaderMixin, SequenceLoaderMixin, OutputWriterMixin, DomainFilterMixin, PredictorMixin
from .._utils import patch_showwarnings

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources  # type: ignore

if typing.TYPE_CHECKING:
    from ...types import TypeClassifier
    from ...model import Cluster, Gene, FeatureTable, GeneTable
    from Bio.SeqIO import SeqRecord


class Predict(TableLoaderMixin, SequenceLoaderMixin, OutputWriterMixin, DomainFilterMixin, PredictorMixin):  # noqa: D101

    summary = "predict gene clusters on contigs that have already been annotated."

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        return f"""
        gecco predict - {cls.summary}

        Usage:
            gecco predict --genome <file> --features <table>... --genes <table> [options]

        Arguments:
            -g <file>, --genome <file>    a genomic file containing one or more
                                          sequences to use as input. Must be in
                                          one of the sequences format supported
                                          by Biopython.
            --features <table>            a feature table obtained by a previous
                                          invocation of ``gecco run``.
            --genes <table>               a gene table obtained by a previous
                                          invocation of ``gecco run``.


        Parameters:
            -f <fmt>, --format <fmt>      the format of the input file, as a
                                          Biopython format string. GECCO is able
                                          to recognize FASTA and GenBank files
                                          automatically if this is not given.
            -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                          multithreading. Use 0 to use all of the
                                          available CPUs. [default: 0]

        Parameters - Output:
            -o <out>, --output-dir <out>  the directory in which to write the
                                          output files. [default: .]
            --antismash-sideload          write an AntiSMASH v6 sideload JSON
                                          file next to the output files.
            --force-tsv                   always write TSV output files even
                                          when they are empty (e.g. because
                                          no genes or no clusters were found).
            --merge-gbk                   output a single file containing
                                          every detected cluster instead of
                                          writing one file per cluster.

        Parameters - Domain Annotation:
            -e <e>, --e-filter <e>        the e-value cutoff for protein domains
                                          to be included. This is not stable
                                          across versions, so consider using
                                          a p-value filter instead.
            -p <p>, --p-filter <p>        the p-value cutoff for protein domains
                                          to be included. [default: 1e-9]

        Parameters - Cluster Detection:
            --no-pad                      disable padding of gene sequences
                                          (used to predict gene clusters in 
                                          contigs smaller than the CRF window 
                                          length).
            -c <N>, --cds <N>             the minimum number of coding sequences a
                                          valid cluster must contain. [default: 3]
            -m <m>, --threshold <m>       the probability threshold for cluster
                                          detection. Default depends on the
                                          post-processing method (0.8 for gecco,
                                          0.6 for antismash).
            --postproc <method>           the method to use for cluster validation
                                          (antismash or gecco). [default: gecco]
            -E <N>, --edge-distance <N>   the minimum number of annotated genes
                                          that must separate a cluster from the
                                          edge. Edge clusters will still be
                                          included if they are longer. A lower
                                          number will increase the number of
                                          false positives on small contigs.
                                          [default: 0]

        Parameters - Debug:
            --model <directory>           the path to an alternative CRF model
                                          to use (obtained with `gecco train`).

        """

    def _check(self) -> None:
        Command._check(self)
        try:
            self.cds = self._check_flag("--cds", int, lambda x: x > 0, hint="positive integer")
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
            if self.args["--threshold"] is None:
                self.threshold = 0.8 if self.args["--postproc"] == "gecco" else 0.6
            else:
                self.threshold = self._check_flag("--threshold", float, lambda x: 0 <= x <= 1, hint="number between 0 and 1")
            self.jobs = self._check_flag("--jobs", int, lambda x: x >= 0, hint="positive or null integer")
            self.postproc = self._check_flag("--postproc", str, lambda x: x in ("gecco", "antismash"), hint="'gecco' or 'antismash'")
            self.edge_distance = self._check_flag("--edge-distance", int, lambda x: x >= 0, hint="positive or null integer")
            self.format = self._check_flag("--format", optional=True)
            self.genome = self._check_flag("--genome")
            self.features: List[str] = self._check_flag("--features")
            self.genes: str = self._check_flag("--genes")
            self.model: Optional[str] = self._check_flag("--model", optional=True)
            self.output_dir: str = self._check_flag("--output-dir")
            self.antismash_sideload = self._check_flag("--antismash-sideload", bool)
            self.force_tsv = self._check_flag("--force-tsv", bool)
            self.no_pad = self._check_flag("--no-pad", bool)
            self.merge_gbk = self._check_flag("--merge-gbk", bool)
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _assign_sources(self, sequences: Iterable["SeqRecord"], genes: List["Gene"]) -> Iterable["Gene"]:
        from ...model import Protein, Strand

        self.info("Extracting", "known sequences from gene list", level=2)
        known_sequences = { gene.source.id for gene in genes }

        self.info("Building", "index of source sequences", level=2)
        sequence_index = {
            record.id:record
            for record in sequences
            if record.id in known_sequences
        }

        try:
            self.info("Assigning", "source sequences to gene objects", level=2)
            for gene in genes:
                gene = gene.with_source(sequence_index[gene.source.id])
                gene_seq = gene.source.seq[gene.start-1:gene.end]
                if gene.strand == Strand.Reverse:
                    gene_seq = gene_seq.reverse_complement()
                gene = gene.with_protein(gene.protein.with_seq(gene_seq.translate()))
                yield gene
        except KeyError as err:
            self.error("Sequence {!r} not found in {!r}", gene.source.id, self.genome)
            raise CommandExit(1) from err

    def _write_sideload_json(self, clusters: List["Cluster"]) -> None:
        # record version and important parameters
        data: Dict[str, Any] = {
            "records": [],
            "tool": {
                "name": "GECCO",
                "version": __version__,
                "description": "Biosynthetic Gene Cluster prediction with Conditional Random Fields.",
                "configuration": {
                    "cds": repr(self.cds),
                    "e-filter": None if self.e_filter is None else repr(self.e_filter),
                    "p-filter": None if self.p_filter is None else repr(self.p_filter),
                    "postproc": repr(self.postproc),
                    "threshold": repr(self.threshold),
                    "edge-distance": repr(self.edge_distance),
                    "no-pad": repr(self.no_pad),
                }
            }
        }
        # record if non-standard HMM or model was used
        if self.model:
            data["tool"]["configuration"]["model"] = self.model
        # create a record per sequence
        for seq_id, seq_clusters in itertools.groupby(clusters, key=operator.attrgetter("source.id")):
            data["records"].append({"name": seq_id, "subregions": []})
            for cluster in seq_clusters:
                probabilities = {
                    f"{key.lower()}_probability":f"{value:.3f}"
                    for key, value in cluster.type_probabilities.items()
                }
                data["records"][-1]["subregions"].append({
                    "start": cluster.start,
                    "end": cluster.end,
                    "label": ";".join(sorted(cluster.type.names)) or "Unknown",
                    "details": {
                        "average_p": f"{cluster.average_probability:.3f}",
                        "max_p": f"{cluster.maximum_probability:.3f}",
                        **probabilities,
                    }
                })
        # write the JSON file to the output folder
        base, _ = os.path.splitext(os.path.basename(self.genome))
        sideload_out = os.path.join(self.output_dir, f"{base}.sideload.json")
        self.info("Writing", "sideload JSON to", repr(sideload_out), level=1)
        with open(sideload_out, "w") as out:
            json.dump(data, out, sort_keys=True, indent=4)

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
            outputs = [f"{base}.clusters.tsv", f"{base}.features.tsv", f"{base}.genes.tsv"]
            if self.antismash_sideload:
                outputs.append(f"{base}.sideload.json")
            self._make_output_directory(outputs)
            # load features
            genes = list(self._load_genes())
            features = self._load_features()
            # load sequences
            sequences = self._load_sequences()
            # label genes
            genes = self._annotate_genes(genes, features)
            genes = list(self._assign_sources(sequences, genes))
            # Sort genes
            self.info("Sorting", "genes by coordinates", level=2)
            genes.sort(key=operator.attrgetter("source.id", "start", "end"))
            for gene in genes:
                gene.protein.domains.sort(key=operator.attrgetter("start", "end"))
            # filter domains by p-value and/or e-value
            genes = self._filter_domains(genes)
            # check if probabilities need to be predicted
            if not all(gene.average_probability is not None for gene in genes):
                genes = self._predict_probabilities(genes)
            self._write_genes_table(genes)
            self._write_feature_table(genes)
            # extract clusters from probability vector
            clusters = self._extract_clusters(genes)
            if clusters:
                self.success("Found", len(clusters), "potential gene clusters", level=1)
            else:
                self.warn("No gene clusters were found")
                if self.force_tsv:
                    self._write_cluster_table(clusters)
                return 0
            # predict types for putative clusters
            classifier = self._load_type_classifier()
            if len(classifier.classes_) > 1:
                clusters = self._predict_types(clusters, classifier)
            # write results
            self.info("Writing", "result files to folder", repr(self.output_dir), level=1)
            self._write_cluster_table(clusters)
            self._write_clusters(clusters, merge=self.merge_gbk)
            if self.antismash_sideload:
                self._write_sideload_json(clusters)
            unit = "cluster" if len(clusters) == 1 else "clusters"
            self.success("Found", len(clusters), "gene", unit, level=0)
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
