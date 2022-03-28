"""Implementation of the ``gecco run`` subcommand.
"""

import contextlib
import errno
import glob
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
from typing import Any, Dict, Union, Optional, List, TextIO, Mapping

import rich.emoji
import rich.progress

from ... import __version__
from ._base import Command, CommandExit, InvalidArgument
from .annotate import Annotate
from .._utils import patch_showwarnings
from ...model import ProductType

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources


class Run(Annotate):  # noqa: D101

    summary = "predict BGC from one or several contigs."

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
        return f"""
        gecco run - {cls.summary}

        Usage:
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
            -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                          multithreading. Use 0 to use all of the
                                          available CPUs. [default: 0]

        Parameters - Output:
            -o <out>, --output-dir <out>  the directory in which to write the
                                          output files. [default: .]
            --antismash-sideload          write an AntiSMASH v6 sideload JSON
                                          file next to the output files.
            --force-clusters-tsv          always write a ``clusters.tsv`` file
                                          even when no clusters were found.

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

        Parameters - Cluster Detection:
            -c <N>, --cds <N>             the minimum number of coding sequences a
                                          valid cluster must contain. [default: 3]
            -m <m>, --threshold <m>       the probability threshold for cluster
                                          detection. Default depends on the
                                          post-processing method (0.3 for gecco,
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
            --hmm <hmm>                   the path to one or more alternative HMM
                                          file to use (in HMMER format).
        """

    def _check(self) -> typing.Optional[int]:
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
                self.threshold = 0.3 if self.args["--postproc"] == "gecco" else 0.6
            else:
                self.threshold = self._check_flag("--threshold", float, lambda x: 0 <= x <= 1, hint="number between 0 and 1")
            self.jobs = self._check_flag("--jobs", int, lambda x: x >= 0, hint="positive or null integer")
            self.postproc = self._check_flag("--postproc", str, lambda x: x in ("gecco", "antismash"), hint="'gecco' or 'antismash'")
            self.edge_distance = self._check_flag("--edge-distance", int, lambda x: x >= 0, hint="positive or null integer")
            self.format = self._check_flag("--format", optional=True)
            self.genome = self._check_flag("--genome")
            self.model = self._check_flag("--model", optional=True)
            self.hmm = self._check_flag("--hmm")
            self.output_dir = self._check_flag("--output-dir")
            self.antismash_sideload = self._check_flag("--antismash-sideload", bool)
            self.force_clusters_tsv = self._check_flag("--force-clusters-tsv", bool)
            self.mask = self._check_flag("--mask", bool)
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _load_model_domains(self) -> typing.Set[str]:
        try:
            if self.model is None:
                self.info("Loading", "feature list from internal model", level=2)
                domains_file = importlib_resources.open_text("gecco.types", "domains.tsv")
            else:
                self.info("Loading", "feature list from", repr(self.model), level=2)
                domains_file = open(os.path.join(self.model, "domains.tsv"))
            with domains_file as f:
                domains = set(filter(None, map(str.strip, f)))
        except FileNotFoundError as err:
            if self.model is not None:
                self.error("Could not find domains list :", repr(self.model))
            raise CommandExit(e.errno) from err
        else:
            self.success("Found", len(domains), "selected features", level=2)
            return domains

    def _predict_probabilities(self, genes):
        from ...crf import ClusterCRF

        if self.model is None:
            self.info("Loading", "embedded CRF pre-trained model", level=1)
        else:
            self.info("Loading", "CRF pre-trained model from", repr(self.model), level=1)
        crf = ClusterCRF.trained(self.model)

        self.info("Predicting", "cluster probabilitites with the CRF model", level=1)
        unit = "genes" if len(genes) > 1 else "gene"
        task = self.progress.add_task("Predicting marginals", total=len(genes), unit=unit, precision="")
        return list(crf.predict_probabilities(
            self.progress.track(genes, task_id=task, total=len(genes)),
            cpus=self.jobs
        ))

    def _extract_clusters(self, genes):
        from ...refine import ClusterRefiner

        self.info("Extracting", "predicted biosynthetic regions", level=1)
        refiner = ClusterRefiner(
            threshold=self.threshold,
            criterion=self.postproc,
            n_cds=self.cds,
            edge_distance=self.edge_distance
        )

        total = len({gene.source.id for gene in genes})
        unit = "contigs" if total > 1 else "contig"
        task = self.progress.add_task("Extracting clusters", total=total, unit=unit, precision="")

        clusters = []
        gene_groups = itertools.groupby(genes, lambda g: g.source.id)
        for _, gene_group in self.progress.track(gene_groups, task_id=task, total=total):
            clusters.extend(refiner.iter_clusters(gene_group))

        return clusters

    def _predict_types(self, clusters):
        from ...model import ProductType
        from ...types import TypeClassifier

        self.info("Predicting", "BGC types", level=1)

        unit = "cluster" if len(clusters) == 1 else "clusters"
        task = self.progress.add_task("Predicting types", total=len(clusters), unit=unit, precision="")

        clusters_new = []
        classifier = TypeClassifier.trained(self.model)
        for cluster in self.progress.track(clusters, task_id=task):
            clusters_new.extend(classifier.predict_types([cluster]))
            if cluster.type != ProductType.Unknown:
                name = "/".join(f"[bold blue]{ty.name}[/]" for ty in cluster.type.unpack())
                prob = "/".join(f"[bold purple]{cluster.type_probabilities[ty]:.0%}[/]" for ty in cluster.type.unpack())
                self.success(f"Predicted type of [bold blue]{cluster.id}[/] as {name} ({prob} confidence)")
            else:
                ty = max(cluster.type_probabilities, key=cluster.type_probabilities.get)
                prob = f"[bold purple]{cluster.type_probabilities[ty]:.0%}[/]"
                name = f"[bold blue]{ty.name}[/]"
                self.warn(f"Couldn't assign type to [bold blue]{cluster.id}[/] (maybe {name}, {prob} confidence)")

        return clusters_new

    def _write_cluster_table(self, clusters):
        from ...model import ClusterTable

        base, _ = os.path.splitext(os.path.basename(self.genome))
        cluster_out = os.path.join(self.output_dir, f"{base}.clusters.tsv")
        self.info("Writing", "cluster table to", repr(cluster_out), level=1)
        with open(cluster_out, "w") as out:
            ClusterTable.from_clusters(clusters).dump(out)

    def _write_clusters(self, clusters):
        from Bio import SeqIO

        for cluster in clusters:
            gbk_out = os.path.join(self.output_dir, f"{cluster.id}.gbk")
            self.info("Writing", f"cluster [bold blue]{cluster.id}[/] to", repr(gbk_out), level=1)
            SeqIO.write(cluster.to_seq_record(), gbk_out, "genbank")

    def _write_sideload_json(self, clusters):
        # record version and important parameters
        data = {
            "records": [],
            "tool": {
                "name": "GECCO",
                "version": __version__,
                "description": "Biosynthetic Gene Cluster prediction with Conditional Random Fields.",
                "configuration": {
                    "cds": repr(self.cds),
                    "e-filter": repr(self.e_filter),
                    "postproc": repr(self.postproc),
                    "threshold": repr(self.threshold),
                    "mask": repr(self.mask),
                    "edge-distance": repr(self.edge_distance),
                }
            }
        }
        # record if non-standard HMM or model was used
        if self.hmm:
            data["tool"]["configuration"]["hmm"] = list(self.hmm)
        if self.model:
            data["tool"]["configuration"]["model"] = self.model
        # create a record per sequence
        for seq_id, seq_clusters in itertools.groupby(clusters, key=lambda cluster: cluster.source.id):
            data["records"].append({"name": seq_id, "subregions": []})
            for cluster in seq_clusters:
                ty = ";".join(sorted(ty.name for ty in cluster.type.unpack())) or "Unknown"
                data["records"][-1]["subregions"].append({
                    "start": cluster.start,
                    "end": cluster.end,
                    "label": ty,
                    "details": {
                        "average_p": f"{cluster.average_probability:.3f}",
                        "max_p": f"{cluster.maximum_probability:.3f}",
                        "alkaloid_probability": f"{cluster.type_probabilities.get(ProductType.Alkaloid, 0.0):.3f}",
                        "polyketide_probability": f"{cluster.type_probabilities.get(ProductType.Polyketide, 0.0):.3f}",
                        "ripp_probability": f"{cluster.type_probabilities.get(ProductType.RiPP, 0.0):.3f}",
                        "saccharide_probability": f"{cluster.type_probabilities.get(ProductType.Saccharide, 0.0):.3f}",
                        "terpene_probability": f"{cluster.type_probabilities.get(ProductType.Terpene, 0.0):.3f}",
                        "nrp_probability": f"{cluster.type_probabilities.get(ProductType.NRP, 0.0):.3f}",
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
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # attempt to create the output directory, checking it doesn't
            # already contain output files (or raise a warning)
            extensions = ["features.tsv", "genes.tsv", "clusters.tsv"]
            if self.antismash_sideload:
                extensions.append("sideload.json")
            self._make_output_directory(extensions)
            # load sequences and extract genes
            sequences = self._load_sequences()
            genes = self._extract_genes(sequences)
            if genes:
                self.success("Found", "a total of", len(genes), "genes", level=1)
            else:
                self.warn("No genes were found")
                return 0
            # use a whitelist for domain annotation, so that we only annotate
            # with features that are useful for the CRF
            whitelist = self._load_model_domains()
            # annotate domains and predict probabilities
            genes = self._annotate_domains(genes, whitelist=whitelist)
            genes = self._predict_probabilities(genes)
            self._write_genes_table(genes)
            self._write_feature_table(genes)
            # extract clusters from probability vector
            clusters = self._extract_clusters(genes)
            if clusters:
                self.success("Found", len(clusters), "potential gene clusters", level=1)
            else:
                self.warn("No gene clusters were found")
                if self.force_clusters_tsv:
                    self._write_cluster_table(clusters)
                return 0
            # predict types for putative clusters
            clusters = self._predict_types(clusters)
            # write results
            self.info("Writing", "result files to folder", repr(self.output_dir), level=1)
            self._write_cluster_table(clusters)
            self._write_clusters(clusters)
            if self.antismash_sideload:
                self._write_sideload_json(clusters)
            unit = "cluster" if len(clusters) == 1 else "clusters"
            self.success("Found", len(clusters), "biosynthetic gene", unit, level=0)
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
