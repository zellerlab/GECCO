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
from typing import Any, Dict, Union, Optional, List, Iterator, TextIO, Mapping, Iterable

from .._utils import in_context, patch_showwarnings, ProgressReader
from ._base import Command, CommandExit, InvalidArgument
from ._mixins import (
    TableLoaderMixin,
    DomainFilterMixin,
    OutputWriterMixin,
    ClusterLoaderMixin,
    CompositionWriterMixin
)

if typing.TYPE_CHECKING:
    from ...crf import ClusterCRF
    from ...model import Cluster, Gene, FeatureTable, ClusterTable


class Train(TableLoaderMixin, DomainFilterMixin, OutputWriterMixin, ClusterLoaderMixin, CompositionWriterMixin):  # noqa: D101

    summary = "train the CRF model on an embedded feature table."

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        return f"""
        gecco train - {cls.summary}

        Usage:
            gecco train --features <table>... [options]

        Arguments:
            -f <data>, --features <table>   a domain annotation table, used to
                                            train the CRF model.
            -c <data>, --clusters <table>   a cluster annotation table, used to
                                            extract the domain composition for
                                            the type classifier.
            -g <file>, --genes <file>       a gene table containing the
                                            coordinates of the genes inside
                                            the training sequence.

        Parameters:
            -j <jobs>, --jobs <jobs>        the number of CPUs to use for
                                            multithreading. Use 0 to use all
                                            the available CPUs. [default: 0]

        Parameters - Output:
            -o <out>, --output-dir <out>    the directory to use for the model
                                            files. [default: model]

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
            -W <N>, --window-size <N>       the size of the sliding window for
                                            CRF predictions. [default: 5]
            --window-step <N>               the step of the sliding window for
                                            CRF predictions. [default: 1]
            --c1 <C1>                       parameter for L1 regularisation.
                                            [default: 0.15]
            --c2 <C2>                       parameter for L2 regularisation.
                                            [default: 0.15]
            --feature-type <type>           at which level features should be
                                            extracted (protein or domain).
                                            [default: protein]
            --select <N>                    fraction of most significant features
                                            to select from the training data.
            --correction <method>           the multiple test correction method
                                            to use when computing significance
                                            with multiple Fisher tests.

        """

    def _check(self) -> None:
        from ...crf.select import _CORRECTION_METHODS

        super()._check()
        try:
            self.feature_type = self._check_flag(
                "--feature-type",
                str,
                lambda x: x in {"protein", "domain"},
                hint="'protein' or 'domain'",
                default="protein",
                optional=True,
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
            self.features: List[str] = self._check_flag("--features", list)
            self.clusters: str = self._check_flag("--clusters", str)
            self.genes: str = self._check_flag("--genes", str)
            self.window_size = self._check_flag(
                "--window-size",
                int,
                lambda x: x > 0,
                hint="positive integer",
            )
            self.window_step = self._check_flag(
                "--window-step",
                int,
                lambda x: x > 0 and x <= self.window_size,
                hint="positive integer smaller than `--window-size`",
            )
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _seed_rng(self) -> None:
        self.info("Seeding", "the random number generator with seed", self.seed, level=2)
        random.seed(self.seed)

    def _fit_model(self, genes: List["Gene"]) -> "ClusterCRF":
        from ...crf import ClusterCRF

        self.info("Creating", f"the CRF in [bold blue]{self.feature_type}[/] mode", level=1)
        self.info("Using", f"provided hyperparameters (C1={self.c1}, C2={self.c2})", level=1)
        if self.select is not None:
            self.info("Selecting", f"features with Fisher Exact Test (threshold={self.select})", level=1)
        self.info("Iterating", f"over features with a sliding window (W={self.window_size}, step={self.window_step})", level=1)
        crf = ClusterCRF(
            self.feature_type,
            algorithm="lbfgs",
            window_size=self.window_size,
            window_step=self.window_step,
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

    # ---

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check arguments and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))  # type: ignore
            # seed RNG
            self._seed_rng()
            # attempt to create the output directory
            outputs = [
                "model.pkl",
                "model.pkl.md5",
                "domains.tsv",
                "types.tsv",
                "compositions.npz"
            ]
            self._make_output_directory(outputs)
            # load features
            genes = list(self._load_genes())
            features = self._load_features()
            genes = self._annotate_genes(genes, features)
            # Sort genes
            self.info("Sorting", "genes by coordinates", level=2)
            genes.sort(key=operator.attrgetter("source.id", "start", "end"))
            for gene in genes:
                gene.protein.domains.sort(key=operator.attrgetter("start", "end"))
            # filter domains by p-value and/or e-value
            genes = self._filter_domains(genes)
            # load clusters and label genes inside clusters
            clusters = self._load_clusters()
            genes = self._label_genes(genes, clusters)
            # fit CRF
            crf = self._fit_model(genes)
            # save model
            self.info("Saving", f"CRF model to {self.output_dir!r}")
            crf.save(self.output_dir)
            self._save_transitions(crf)
            self._save_weights(crf)
            # extract the domains names
            self.info("Finding", "the array of possible protein domains", level=2)
            if crf.significant_features is not None:
                all_possible = sorted(crf.significant_features)
            else:
                all_possible = sorted({d.name for g in genes for d in g.protein.domains})
            # compute domain compositions
            self._save_domain_compositions(all_possible, self._extract_clusters(genes, clusters))
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
