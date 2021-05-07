"""Implementation of the ``gecco cv`` subcommand.
"""

import contextlib
import copy
import functools
import itertools
import os
import operator
import multiprocessing
import random
import typing
from typing import Any, Dict, Union, Optional, List, TextIO, Mapping

from ._base import Command, CommandExit, InvalidArgument
from .._utils import patch_showwarnings


class Cv(Command):  # noqa: D101

    summary = "perform cross validation on a training set."

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
        return f"""
        gecco cv  - {cls.summary}

        Usage:
            gecco cv kfold -f <table> [-c <data>] [options]
            gecco cv loto  -f <table>  -c <data>  [options]

        Arguments:
            -f <data>, --features <table>   a domain annotation table, used to
                                            labeled as BGCs and non-BGCs.
            -c <data>, --clusters <table>   a cluster annotation table, used
                                            to stratify clusters by type in
                                            LOTO mode.

        Parameters:
            -o <out>, --output <out>        the name of the output file where
                                            the cross-validation table will be
                                            written. [default: cv.tsv]
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
            --overlap <N>                   how much overlap to consider if
                                            features overlap. [default: 2]
            --splits <N>                    number of folds for cross-validation
                                            (if running `kfold`). [default: 10]
            --select <N>                    fraction of most significant features
                                            to select from the training data.
            --shuffle                       enable shuffling of stratified rows.

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
            self.overlap = self._check_flag(
                "--overlap",
                int,
                lambda x: x > 0,
                hint="positive integer",
            )
            self.c1 = self._check_flag("--c1", float, hint="real number")
            self.c2 = self._check_flag("--c2", float, hint="real number")
            self.splits = self._check_flag("--splits", int, lambda x: x>1, hint="integer greater than 1")
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
            self.clusters = self._check_flag("--clusters")
            self.loto = self.args["loto"]
            self.output = self.args["--output"]
        except InvalidArgument:
            raise CommandExit(1)

    # --

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

    def _group_genes(self, genes):
        self.info("Grouping", "genes by source sequence")
        groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
        seqs = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]
        if self.args["--shuffle"]:
            self.info("Shuffling training data sequences")
            random.shuffle(seqs)
        return seqs

    def _loto_splits(self, seqs):
        from ...crf.cv import LeaveOneGroupOut
        from ...model import ClusterTable, ProductType

        self.info("Loading", "the clusters table")
        with open(self.clusters) as in_:
            table = ClusterTable.load(in_)
            index = { row.sequence_id: row.type for row in table }
            if len(index) != len(table):
                raise ValueError("Training data contains several clusters per sequence")

        self.info("Grouping", "sequences by cluster types")
        groups = []
        for cluster in seqs:
            ty = next((index[g.source.id] for g in cluster if g.source.id in index), None)
            if ty is None:
                seq_id = next(gene.source.id for gene in cluster)
                self.logger.warning("Could not find type of cluster in {!r}", seq_id)
                ty = ProductType.Unknown
            groups.append(ty.unpack())

        return list(LeaveOneGroupOut().split(seqs, groups=groups))

    def _kfold_splits(self, seqs):
        import sklearn.model_selection
        return list(sklearn.model_selection.KFold(self.splits).split(seqs))

    def _get_train_data(self, train_indices, seqs):
        # extract train data
        return [gene for i in train_indices for gene in seqs[i]]

    def _get_test_data(self, test_indices, seqs):
        # make a clean copy of the test data without gene probabilities
        test_data = [copy.deepcopy(gene) for i in test_indices for gene in seqs[i]]
        for gene in test_data:
            gene.protein.domains = [d.with_probability(None) for d in gene.protein.domains]
        return test_data

    def _fit_predict(self, train_data, test_data):
        from ...crf import ClusterCRF

        # fit and predict the CRF for the current fold
        crf = ClusterCRF(
            self.feature_type,
            algorithm="lbfgs",
            overlap=self.overlap,
            c1=self.c1,
            c2=self.c2,
        )
        crf.fit(train_data, cpus=self.jobs, select=self.select)
        return crf.predict_probabilities(test_data, cpus=self.jobs)

    def _write_fold(self, fold, genes, append=False):
        from ...model import FeatureTable

        frame = FeatureTable.from_genes(genes).to_dataframe()
        with open(self.output, "a" if append else "w") as out:
            frame.assign(fold=fold).to_csv(out, header=not append, sep="\t", index=False)

    # --

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check arguments and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # load features
            features = self._load_features()
            self.success("Loaded", len(features), "feature annotations")
            genes = self._convert_to_genes(features)
            self.success("Recoverd", len(genes), "genes from the feature annotations")
            seqs = self._group_genes(genes)
            self.success("Grouped", "genes into", len(seqs), "sequences")
            # split CV folds
            if self.loto:
                splits = self._loto_splits(seqs)
            else:
                splits = self._kfold_splits(seqs)
            # run CV
            unit = "fold" if len(splits) == 1 else "folds"
            task = self.progress.add_task(description="Cross-Validation", total=len(splits), unit=unit)
            self.info("Performing cross-validation")
            for i, (train_indices, test_indices) in enumerate(self.progress.track(splits, task_id=task)):
                train_data = self._get_train_data(train_indices, seqs)
                test_data = self._get_test_data(test_indices, seqs)
                new_genes = self._fit_predict(train_data, test_data)
                self._write_fold(i+1, new_genes, append=i>0)
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
