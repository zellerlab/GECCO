"""Implementation of the ``gecco train`` subcommand.
"""

import contextlib
import csv
import hashlib
import io
import os
import operator
import pickle
import signal
import typing
from typing import Any, Dict, Union, Optional, List, TextIO, Mapping

from ._base import Command, CommandExit, InvalidArgument
from .._utils import in_context, patch_showwarnings


class Train(Command):  # noqa: D101

    summary = "train the CRF model on an embedded feature table."

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
        return f"""
        gecco train - {cls.summary}

        Usage:
            gecco train --features <table> --clusters <table> [options]

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
                                            be included [default: 1e-5]

        Parameters - Training:
            --c1 <C1>                       parameter for L1 regularisation.
                                            [default: 0.15]
            --c2 <C2>                       parameter for L2 regularisation.
                                            [default: 0.15]
            --feature-type <type>           how features should be extracted
                                            (single, overlap, or group).
                                            [default: group]
            --truncate <N>                  the maximum number of rows to use from
                                            the training set.
            --overlap <N>                   how much overlap to consider if
                                            features overlap. [default: 2]
            --no-shuffle                    disable shuffling of the data before
                                            fitting the model.
            --select <N>                    fraction of most significant features
                                            to select from the training data.

        """

    def _check(self) -> typing.Optional[int]:
        super()._check()
        try:
            self.feature_type = self._check_flag(
                "--feature-type",
                str,
                lambda x: x in {"single", "overlap", "group"},
                hint="'single', 'overlap' or 'group'"
            )
            self.truncate = self._check_flag(
                "--truncate",
                lambda x: x if x is None else int(x),
                lambda x: x is None or x > 0,
                hint="positive integer"
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
                lambda x: 0 <= x <= 1,
                hint="real number between 0 and 1"
            )
            self.select = self._check_flag(
                "--select",
                lambda x: x if x is None else float(x),
                lambda x: x is None or 0 <= x <= 1,
                hint="real number between 0 and 1"
            )
            self.jobs = self._check_flag(
                "--jobs",
                int,
                lambda x: x >= 0,
                hint="positive or null integer"
            )
            self.features = self._check_flag("--features")
            self.no_shuffle = self._check_flag("--no-shuffle", bool)
            self.output_dir = self._check_flag("--output-dir", str)
            self.features = self._check_flag("--features", str)
            self.clusters = self._check_flag("--clusters", str)
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _make_output_directory(self) -> None:
        # Make output directory
        self.info("Using", "output folder", repr(self.output_dir), level=1)
        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except OSError as err:
            self.error("Could not create output directory: {}", err)
            raise CommandExit(e.errno) from err
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

    def _load_features(self):
        from ...model import FeatureTable

        self.info("Loading", "features table from file", repr(self.features))
        with open(self.features) as in_:
            return FeatureTable.load(in_)

    def _convert_to_genes(self, features):
        self.info("Converting", "features to genes")

        gene_count = len(set(features.protein_id))
        unit = "gene" if gene_count == 1 else "genes"
        task = self.progress.add_task("Feature conversion", total=gene_count, unit=unit)

        genes = list(self.progress.track(
            features.to_genes(),
            total=gene_count,
            task_id=task
        ))

        self.info("Sorting", "genes by genomic coordinates")
        genes.sort(key=operator.attrgetter("source.id", "start", "end"))
        self.info("Sorting", "domains by protein coordinates")
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))
        return genes

    def _fit_model(self, genes):
        from ...crf import ClusterCRF

        self.info("Creating" f"the CRF in {self.feature_type} mode", level=2)
        self.info("Using" f"hyperparameters C1={self.c1}, C2={self.c2}", level=2)
        crf = ClusterCRF(
            self.feature_type,
            algorithm="lbfgs",
            overlap=self.overlap,
            c1=self.c1,
            c2=self.c2,
        )
        self.info("Fitting", "the CRF model to the training data")
        crf.fit(genes, select=self.select, shuffle=not self.no_shuffle, cpus=self.jobs)
        return crf

    def _save_model(self, crf):
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

    def _save_transitions(self, crf):
        self.info("Writing", "CRF transitions weights")
        with open(os.path.join(self.output_dir, "model.trans.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["from", "to", "weight"])
            for labels, weight in crf.model.transition_features_.items():
                writer.writerow([*labels, weight])

    def _save_weights(self, crf):
        self.info("Writing", "state weights")
        with open(os.path.join(self.output_dir, "model.state.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["attr", "label", "weight"])
            for attrs, weight in crf.model.state_features_.items():
                writer.writerow([*attrs, weight])

    def _load_clusters(self, genes):
        from ...model import ClusterTable, Cluster

        self.info("Loading", "clusters table from file", repr(self.clusters))
        index = { g.id: g for g in genes }
        with open(self.clusters) as f:
            clusters = [
                Cluster(c.bgc_id, [index[p] for p in c.proteins if p in index ], c.type)
                for c in ClusterTable.load(f)
            ]
        clusters.sort(key=operator.attrgetter("id"))
        return clusters

    def _save_domain_compositions(self, crf, clusters):
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
        comp = numpy.array([c.domain_composition(all_possible) for c in clusters])
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
            # attempt to create the output directory
            self._make_output_directory()
            # load features
            features = self._load_features()
            genes = self._convert_to_genes(features)
            del features
            # fit CRF
            crf = self._fit_model(genes)
            # save model
            self._save_model(crf)
            self._save_transitions(crf)
            self._save_weights(crf)
            # load clusters
            clusters = self._load_clusters(genes)
            self._save_domain_compositions(crf, clusters)
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
