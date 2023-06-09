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
from ._mixins import SequenceLoaderMixin, OutputWriterMixin, PredictorMixin
from .._utils import patch_showwarnings

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files # type: ignore

if typing.TYPE_CHECKING:
    from ...types import TypeClassifier
    from ...model import Cluster, Gene


class Run(Annotate, SequenceLoaderMixin, OutputWriterMixin, PredictorMixin):  # noqa: D101

    summary = "predict gene clusters from one or several contigs."

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        return f"""
        gecco run - {cls.summary}

        Usage:
            gecco run --genome <file> [--hmm <hmm>]... [--hmm-x <hmm>]... [options]

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
            --force-tsv                   always write TSV output files even
                                          when they are empty (e.g. because
                                          no genes or no clusters were found).
            --merge-gbk                   output a single file containing
                                          every detected cluster instead of
                                          writing one file per cluster.

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

        Parameters - Cluster Detection:
            --model <directory>           the path to an alternative CRF model
                                          to use (obtained with `gecco train`).
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

        """

    def _check(self) -> None:
        _BITSCORE_CUTOFFS = {"gathering", "noise", "trusted"}
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
            self.model: Optional[str] = self._check_flag("--model", optional=True)
            self.hmm = self._check_flag("--hmm")
            self.output_dir = self._check_flag("--output-dir")
            self.antismash_sideload = self._check_flag("--antismash-sideload", bool)
            self.force_tsv = self._check_flag("--force-tsv", bool)
            self.mask = self._check_flag("--mask", bool)
            self.cds_feature = self._check_flag("--cds-feature", optional=True)
            self.locus_tag = self._check_flag("--locus-tag")
            self.no_pad = self._check_flag("--no-pad", bool)
            self.merge_gbk = self._check_flag("--merge-gbk", bool)
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

    def _load_model_domains(self) -> typing.Set[str]:
        try:
            if self.model is None:
                self.info("Loading", "feature list from internal model", level=2)
                domains_file = files("gecco.types").joinpath("domains.tsv").open()
            else:
                self.info("Loading", "feature list from", repr(self.model), level=2)
                domains_file = open(os.path.join(self.model, "domains.tsv"))
            with domains_file as f:
                domains = set(filter(None, map(str.strip, f)))
        except FileNotFoundError as err:
            if self.model is not None:
                self.error("Could not find domains list :", repr(self.model))
            raise CommandExit(err.errno) from err
        else:
            self.success("Found", len(domains), "selected features", level=2)
            return domains

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
                    "e-filter": repr(self.e_filter),
                    "p-filter": repr(self.p_filter),
                    "bit-cutoffs": repr(self.bit_cutoffs),
                    "postproc": repr(self.postproc),
                    "threshold": repr(self.threshold),
                    "mask": repr(self.mask),
                    "edge-distance": repr(self.edge_distance),
                    "no-pad": repr(self.no_pad),
                }
            }
        }
        # record if non-standard HMM or model was used
        if self.hmm:
            data["tool"]["configuration"]["hmm"] = list(self.hmm)
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
            # load sequences and extract genes
            sequences = list(self._load_sequences())
            genes = self._extract_genes(sequences)
            if genes:
                self.success("Found", "a total of", len(genes), "genes", level=1)
            else:
                if self.force_tsv:
                    self._write_genes_table(genes)
                    self._write_feature_table([])
                    self._write_cluster_table([])
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
