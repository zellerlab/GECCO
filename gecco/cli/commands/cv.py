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
import signal
import typing
from typing import Any, Dict, Union, Optional, Iterable, List, TextIO, Mapping, Tuple

import docopt

from ..._meta import zopen
from .._utils import patch_showwarnings
from ._base import Command, CommandExit, InvalidArgument
from .annotate import Annotate
from .train import Train

if typing.TYPE_CHECKING:
    import numpy
    from numpy.typing import NDArray
    from Bio.SeqRecord import SeqRecord
    from ...hmmer import HMM
    from ...model import Gene
    from ...orf import ORFFinder


class Cv(Train):  # noqa: D101

    summary = "perform cross validation on a training set."

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        return f"""
        gecco cv  - {cls.summary}

        Usage:
            gecco cv kfold --features <table>... --clusters <table> [options]
            gecco cv loto  --features <table>... --clusters <table> [options]

        Arguments:
            -f <data>, --features <table>   a domain annotation table, used
                                            to train the CRF.
            -c <data>, --clusters <table>   a cluster annotation table, used
                                            to stratify clusters by type in
                                            LOTO mode.
            -g <file>, --genes <file>       a gene table containing the
                                            coordinates of the genes inside
                                            the training sequence.

        Parameters:
            -o <out>, --output <out>        the name of the output file where
                                            the cross-validation table will be
                                            written. [default: cv.tsv]
            -j <jobs>, --jobs <jobs>        the number of CPUs to use for
                                            multithreading. Use 0 to use all
                                            the available CPUs. [default: 0]

        Parameters - Domain Annotation:
            -e <e>, --e-filter <e>          the e-value cutoff for domains to
                                            be included.
            -p <p>, --p-filter <p>          the p-value cutoff for domains
                                            to be included. [default: 1e-9]

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

        Parameters - Cross-validation:
            --splits <N>                    number of folds for cross-validation
                                            (if running `kfold`). [default: 10]

        """

    def _check(self) -> None:
        if not isinstance(self.args, docopt.DocoptExit):
            self.args["--output-dir"] = "."
        super()._check()
        try:
            self.output = self._check_flag("--output", str)
            self.splits = self._check_flag(
                "--splits",
                int,
                lambda x: x > 0,
                hint="positive integer"
            )
            self.loto = self.args["loto"]
        except InvalidArgument:
            raise CommandExit(1)

    # --

    def _group_genes(self, genes: List["Gene"]) -> List[List["Gene"]]:
        self.info("Grouping", "genes by source sequence")
        groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
        seqs = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]
        if not self.no_shuffle:
            self.info("Shuffling", "training data sequences")
            random.shuffle(seqs)
        return seqs

    def _loto_splits(self, seqs: List[List["Gene"]]) -> List[Tuple["NDArray[numpy.bool_]", "NDArray[numpy.bool_]"]]:
        from ...crf.cv import LeaveOneGroupOut
        from ...model import ClusterTable, ClusterType

        self.info("Loading", "the clusters table")
        with zopen(self.clusters) as in_:
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
                self.warn("Failed", f"to find type of cluster in {seq_id!r}")
                ty = ClusterType()
            groups.append(ty.unpack())

        return list(LeaveOneGroupOut().split(seqs, groups=groups))   # type: ignore

    def _kfold_splits(self, seqs: List[List["Gene"]]) -> List[Tuple["NDArray[numpy.bool_]", "NDArray[numpy.bool_]"]]:
        import sklearn.model_selection
        return list(sklearn.model_selection.KFold(self.splits).split(seqs))

    @staticmethod
    def _get_train_data(train_indices: Iterable[int], seqs: List[List["Gene"]]) -> List["Gene"]:
        # extract train data
        return [gene for i in train_indices for gene in seqs[i]]

    @staticmethod
    def _get_test_data(test_indices: Iterable[int], seqs: List[List["Gene"]]) -> List["Gene"]:
        # make a clean copy of the test data without gene probabilities
        return [
            gene.with_protein(gene.protein.with_domains(
                d.with_probability(None) for d in gene.protein.domains
            ))
            for i in test_indices
            for gene in seqs[i]
        ]

    def _fit_predict(self, train_data: List["Gene"], test_data: List["Gene"]) -> List["Gene"]:
        from ...crf import ClusterCRF

        crf = self._fit_model(train_data)
        return crf.predict_probabilities(test_data)

    def _write_fold(
        self, 
        fold: int, 
        truth: List["Genes"], 
        predicted: List["Gene"], 
        append: bool = False,
    ) -> None:
        import polars
        from ...model import GeneTable

        frame = (
            GeneTable
                .from_genes(predicted)
                .data
                .with_columns(polars.lit(fold).alias("fold"))
                .with_columns(
                    GeneTable
                        .from_genes(truth)
                        .data
                        .select(is_cluster=polars.col("average_p") > 0.5)
                )
        )
        with open(self.output, "ab" if append else "wb") as out:
            frame.write_csv(out, has_header=not append, separator="\t")

    def _report_fold(
        self,
        fold: typing.Optional[int],
        truth: List["Genes"],
        predicted: List["Genes"],
    ) -> None:
        from sklearn.metrics import average_precision_score, roc_auc_score

        probas = [ gene.average_probability for gene in predicted ]
        labels = [ gene.average_probability > 0.5 for gene in truth ]

        aupr = average_precision_score(labels, probas)
        auroc = roc_auc_score(labels, probas)
        if fold:
            self.info(f"Finished training fold {fold} (AUROC={auroc:.3f}, AUPR={aupr:.3f})")
        else:
            self.info(f"Finished cross validation (AUROC={auroc:.3f}, AUPR={aupr:.3f})")

    # --

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check arguments and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))   # type: ignore
            # seed RNG
            self._seed_rng()
            # load features
            genes = list(self._load_genes())
            features = self._load_features()
            genes = self._annotate_genes(genes, features)
            # load clusters and label genes inside clusters
            clusters = self._load_clusters()
            genes = self._label_genes(genes, clusters)
            seqs = self._group_genes(genes)
            self.success("Grouped", "genes into", len(seqs), "sequences")
            # split CV folds
            if self.loto:
                splits = self._loto_splits(seqs)
            else:
                splits = self._kfold_splits(seqs)
            # run CV
            unit = "fold" if len(splits) == 1 else "folds"
            task = self.progress.add_task(description="Cross-validating", total=len(splits), unit=unit, precision="")
            self.info("Performing cross-validation")
            predicted = []
            for i, (train_indices, test_indices) in enumerate(self.progress.track(splits, task_id=task)):
                train_data = self._get_train_data(train_indices, seqs)
                test_data = self._get_test_data(test_indices, seqs)
                old_genes = self._get_train_data(test_indices, seqs)
                new_genes = self._fit_predict(train_data, test_data)
                self._write_fold(i+1, old_genes, new_genes, append=i>0)
                self._report_fold(i+1, old_genes, new_genes)
                predicted.extend(new_genes)
            self._report_fold(None, genes, predicted)
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
