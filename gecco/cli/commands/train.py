"""Implementation of the ``gecco train`` subcommand.
"""

import collections
import contextlib
import csv
import hashlib
import io
import itertools
import os
import operator
import pickle
import random
import signal
import typing
from typing import Any, Dict, Union, Optional, List, TextIO, Mapping

from .._utils import in_context, patch_showwarnings, ProgressReader
from ._base import Command, CommandExit, InvalidArgument
from .annotate import Annotate

if typing.TYPE_CHECKING:
    from ...crf import ClusterCRF
    from ...model import Cluster, Gene, FeatureTable, ClusterTable


class Train(Command):  # noqa: D101

    summary = "train the CRF model on an embedded feature table."

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        return f"""
        gecco train - {cls.summary}

        Usage:
            gecco train --features <table>... --clusters <table> [options]

        Arguments:
            -f <data>, --features <table>   a domain annotation table, used to
                                            train the CRF model.
            -c <data>, --clusters <table>   a cluster annotation table, used to
                                            extract the domain composition for
                                            the type classifier.

        Parameters:
            -o <out>, --output-dir <out>    the directory to use for the model
                                            files. [default: model]
            -j <jobs>, --jobs <jobs>        the number of CPUs to use for
                                            multithreading. Use 0 to use all
                                            the available CPUs. [default: 0]

        Parameters - Domain Annotation:
            -e <e>, --e-filter <e>          the e-value cutoff for domains to
                                            be included.
            -p <p>, --p-filter <p>          the p-value cutoff for domains to
                                            be included. [default: 1e-9]

        Parameters - Training Data:
            --no-shuffle                    disable shuffling of the data
                                            before fitting the model.
            --seed <N>                      the seed to initialize the RNG
                                            with for shuffling operations.
                                            [default: 42]

        Parameters - Training:
            --c1 <C1>                       parameter for L1 regularisation.
                                            [default: 0.15]
            --c2 <C2>                       parameter for L2 regularisation.
                                            [default: 0.15]
            --feature-type <type>           how features should be extracted
                                            (single, overlap, or group).
                                            [default: group]
            --overlap <N>                   how much overlap to consider if
                                            features overlap. [default: 2]
            --select <N>                    fraction of most significant features
                                            to select from the training data.
            --correction <method>           the multiple test correction method
                                            to use when computing significance
                                            with multiple Fisher tests.

        """

    def _check(self) -> typing.Optional[int]:
        from ...crf.select import _CORRECTION_METHODS

        super()._check()
        try:
            self.feature_type = self._check_flag(
                "--feature-type",
                str,
                lambda x: x in {"single", "overlap", "group"},
                hint="'single', 'overlap' or 'group'"
            )
            self.overlap = self._check_flag(
                "--overlap",
                int,
                lambda x: x > 0,
                hint="positive integer",
            )
            self.c1 = self._check_flag("--c1", float, hint="real number")
            self.c2 = self._check_flag("--c2", float, hint="real number")
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
            self.select = self._check_flag(
                "--select",
                float,
                lambda x: 0 <= x <= 1,
                hint="real number between 0 and 1",
                optional=True,
            )
            self.correction = self._check_flag(
                "--correction",
                str,
                lambda m: m in _CORRECTION_METHODS,
                hint="one of {}".format(", ".join(sorted(_CORRECTION_METHODS))),
                optional=True
            )
            self.jobs = self._check_flag(
                "--jobs",
                int,
                lambda x: x >= 0,
                hint="positive or null integer"
            )
            self.no_shuffle = self._check_flag("--no-shuffle", bool)
            self.seed = self._check_flag("--seed", int)
            self.output_dir = self._check_flag("--output-dir", str)
            self.features = self._check_flag("--features", list)
            self.clusters = self._check_flag("--clusters", str)
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _seed_rng(self):
        self.info("Seeding", "the random number generator with seed", self.seed, level=2)
        random.seed(self.seed)

    def _make_output_directory(self) -> None:
        # Make output directory
        self.info("Using", "output folder", repr(self.output_dir), level=1)
        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except OSError as err:
            self.error("Could not create output directory: {}", err)
            raise CommandExit(err.errno) from err
        # Check if output files already exist
        files = [
            "model.pkl",
            "model.pkl.md5",
            "domains.tsv",
            "types.tsv",
            "compositions.npz"
        ]
        for f in files:
            if os.path.isfile(os.path.join(self.output_dir, f)):
                self.warn("Output folder contains files that will be overwritten")
                break

    def _load_features(self) -> "FeatureTable":
        from ...model import FeatureTable

        features = FeatureTable()
        for filename in self.features:
            try:
                # get filesize and unit
                input_size = os.stat(filename).st_size
                total, scale, unit = ProgressReader.scale_size(input_size)
                task = self.progress.add_task("Loading features", total=total, unit=unit, precision=".1f")
                # load features
                self.info("Loading", "features table from file", repr(filename))
                with ProgressReader(open(filename, "rb"), self.progress, task, scale) as in_:
                    features += FeatureTable.load(io.TextIOWrapper(in_))
            except FileNotFoundError as err:
                self.error("Could not find feature file:", repr(filename))
                raise CommandExit(err.errno) from err

        self.success("Loaded", "a total of", len(features), "features", level=1)
        return features

    def _convert_to_genes(self, features: "FeatureTable") -> List["Gene"]:
        self.info("Converting", "features to genes")

        gene_count = len(set(features.protein_id))
        unit = "gene" if gene_count == 1 else "genes"
        task = self.progress.add_task("Converting features", total=gene_count, unit=unit, precision="")

        genes = list(self.progress.track(
            features.to_genes(),
            total=gene_count,
            task_id=task
        ))

        # filter domains out
        Annotate._filter_domains(self, genes)

        self.info("Sorting", "genes by genomic coordinates")
        genes.sort(key=operator.attrgetter("source.id", "start", "end"))
        self.info("Sorting", "domains by protein coordinates")
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))
        return genes

    def _fit_model(self, genes: List["Gene"]) -> "ClusterCRF":
        from ...crf import ClusterCRF

        self.info("Creating", f"the CRF in [bold blue]{self.feature_type}[/] mode", level=1)
        self.info("Using", f"hyperparameters C1={self.c1}, C2={self.c2}", level=1)
        if self.select is not None:
            self.info("Using", f"Fisher Exact Test significance threshold of {self.select}", level=1)
        crf = ClusterCRF(
            self.feature_type,
            algorithm="lbfgs",
            overlap=self.overlap,
            c1=self.c1,
            c2=self.c2,
        )
        self.info("Fitting", "the CRF model to the training data")
        crf.fit(
            genes,
            select=self.select,
            shuffle=not self.no_shuffle,
            correction_method=self.correction,
            cpus=self.jobs
        )
        return crf

    def _save_model(self, crf: "ClusterCRF") -> None:
        model_out = os.path.join(self.output_dir, "model.pkl")
        self.info("Pickling", "the model to", repr(model_out))
        with open(model_out, "wb") as out:
            pickle.dump(crf, out, protocol=4)

        self.info("Computing", "pickled model checksum", level=2)
        hasher = hashlib.md5()
        with open(model_out, "rb") as out:
            for chunk in iter(lambda: out.read(io.DEFAULT_BUFFER_SIZE), b""):
                hasher.update(chunk)

        self.info("Writing", "pickled model checksum to", repr(f"{model_out}.md5"), level=2)
        with open(f"{model_out}.md5", "w") as out_hash:
            out_hash.write(hasher.hexdigest())

    def _save_transitions(self, crf: "ClusterCRF") -> None:
        self.info("Writing", "CRF transitions weights")
        with open(os.path.join(self.output_dir, "model.trans.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["from", "to", "weight"])
            for labels, weight in crf.model.transition_features_.items():
                writer.writerow([*labels, weight])

    def _save_weights(self, crf: "ClusterCRF") -> None:
        self.info("Writing", "state weights")
        with open(os.path.join(self.output_dir, "model.state.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["attr", "label", "weight"])
            for attrs, weight in crf.model.state_features_.items():
                writer.writerow([*attrs, weight])

    def _load_clusters(self) -> "ClusterTable":
        from ...model import ClusterTable

        try:
            # get filesize and unit
            input_size = os.stat(self.clusters).st_size
            total, scale, unit = ProgressReader.scale_size(input_size)
            task = self.progress.add_task("Loading clusters", total=total, unit=unit, precision=".1f")
            # load clusters
            self.info("Loading", "clusters table from file", repr(self.clusters))
            with ProgressReader(open(self.clusters, "rb"), self.progress, task, scale) as in_:
                return ClusterTable.load(io.TextIOWrapper(in_))
        except FileNotFoundError as err:
            self.error("Could not find clusters file:", repr(self.clusters))
            raise CommandExit(err.errno) from err

    def _label_genes(self, genes: List["Gene"], clusters: "ClusterTable") -> List["Gene"]:
        cluster_by_seq = collections.defaultdict(list)
        for cluster_row in clusters:
            cluster_by_seq[cluster_row.sequence_id].append(cluster_row)

        gene_count = len(genes)
        unit = "gene" if gene_count == 1 else "genes"
        task = self.progress.add_task("Labelling genes", total=gene_count, unit=unit, precision="")

        self.info("Labelling", "genes belonging to clusters")
        labelled_genes = []
        for seq_id, seq_genes in itertools.groupby(genes, key=lambda g: g.source.id):
            for gene in seq_genes:
                if any(
                    cluster_row.start <= gene.start and gene.end <= cluster_row.end
                    for cluster_row in cluster_by_seq[seq_id]
                ):
                    gene.protein.domains = [d.with_probability(1) for d in gene.protein.domains]
                else:
                    gene.protein.domains = [d.with_probability(0) for d in gene.protein.domains]
                labelled_genes.append(gene)
                self.progress.update(task_id=task, advance=1)

        return labelled_genes

    def _extract_clusters(self, genes: List["Gene"], clusters: "ClusterTable") -> List["Cluster"]:
        from ...model import Cluster

        cluster_by_seq = collections.defaultdict(list)
        for cluster_row in clusters:
            cluster_by_seq[cluster_row.sequence_id].append(cluster_row)

        self.info("Extracting", "genes belonging to clusters")
        genes_by_cluster = collections.defaultdict(list)
        for seq_id, seq_genes in itertools.groupby(genes, key=lambda g: g.source.id):
            for gene in seq_genes:
                for cluster_row in cluster_by_seq[seq_id]:
                    if cluster_row.start <= gene.start and gene.end <= cluster_row.end:
                        genes_by_cluster[cluster_row.bgc_id].append(gene)

        return [
            Cluster(cluster_row.bgc_id, genes_by_cluster[cluster_row.bgc_id], cluster_row.type)
            for cluster_row in sorted(clusters, key=lambda row: row.bgc_id)
            if genes_by_cluster[cluster_row.bgc_id]
        ]

    def _save_domain_compositions(self, crf: "ClusterCRF", clusters: List["Cluster"]):
        import numpy
        import scipy.sparse

        self.info("Finding", "the array of possible protein domains", level=2)
        if crf.significant_features is not None:
            all_possible = sorted(crf.significant_features)
        else:
            all_possible = sorted({d.name for c in clusters for g in c.genes for d in g.protein.domains})

        self.info("Saving", "training matrix labels for BGC type classifier")
        with open(os.path.join(self.output_dir, "domains.tsv"), "w") as out:
            out.writelines(f"{domain}\n" for domain in all_possible)
        with open(os.path.join(self.output_dir, "types.tsv"), "w") as out:
            writer = csv.writer(out, dialect="excel-tab")
            for cluster in clusters:
                types =  ";".join(ty.name for ty in cluster.type.unpack())
                writer.writerow([cluster.id, types])

        self.info("Building", "new domain composition matrix")
        comp = numpy.array([
            c.domain_composition(all_possible)
            for c in clusters
        ])

        comp_out = os.path.join(self.output_dir, "compositions.npz")
        self.info("Saving", "new domain composition matrix to file", repr(comp_out))
        scipy.sparse.save_npz(comp_out, scipy.sparse.coo_matrix(comp))

    # ---

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check arguments and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # seed RNG
            self._seed_rng()
            # attempt to create the output directory
            self._make_output_directory()
            # load features
            features = self._load_features()
            genes = self._convert_to_genes(features)
            del features
            # load clusters and label genes inside clusters
            clusters = self._load_clusters()
            genes = self._label_genes(genes, clusters)
            # fit CRF
            crf = self._fit_model(genes)
            # save model
            self._save_model(crf)
            self._save_transitions(crf)
            self._save_weights(crf)
            # compute domain compositions
            self._save_domain_compositions(crf, self._extract_clusters(genes, clusters))
            self.success("Finished", "training new CRF model", level=0)
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
