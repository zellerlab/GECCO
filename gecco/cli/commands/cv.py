"""Implementation of the ``gecco cv`` subcommand.
"""

import functools
import itertools
import os
import operator
import multiprocessing
import random
import typing
from typing import List

import tqdm
import sklearn.model_selection

from ._base import Command
from ...model import FeatureTable
from ...crf import ClusterCRF
from ...crf.cv import LeaveOneGroupOut


class Cv(Command):  # noqa: D101

    summary = "perform cross validation on a training set."
    doc = f"""
    gecco cv  - {summary}

    Usage:
        gecco cv (-h | --help)
        gecco cv kfold -f <table> [-c <data>] [options]
        gecco cv loto  -f <table>  -c <data>  [options]

    Arguments:
        -f <data>, --features <table>   a domain annotation table, used to
                                        labeled as BGCs and non-BGCs.
        -c <data>, --clusters <table>   a cluster annotation table, use to
                                        stratify clusters by type in LOTO
                                        mode.

    Parameters:
        -o <out>, --output <out>        the name of the output cross-validation
                                        table. [default: cv.tsv]
        -j <jobs>, --jobs <jobs>        the number of CPUs to use for
                                        multithreading. Use 0 to use all of the
                                        available CPUs. [default: 0]

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
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check the inputs exist
        for input_ in filter(None, (self.args["--features"], self.args["--clusters"])):
            if not os.path.exists(input_):
                self.logger.error("could not locate input file: {!r}", input_)
                return 1

        # Check the `--feature-type`
        type_ = self.args["--feature-type"]
        if type_ not in {"single", "overlap", "group"}:
            self.logger.error("Invalid value for `--feature-type`: {}", type_)
            return 1

        # Check value of numeric arguments
        self.args["--overlap"] = int(self.args["--overlap"])
        self.args["--c1"] = float(self.args["--c1"])
        self.args["--c2"] = float(self.args["--c2"])
        self.args["--splits"] = int(self.args["--splits"])
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1
        if self.args["--select"] is not None:
            self.args["--select"] = float(self.args["--select"])

        # Check the `--jobs`flag
        self.args["--jobs"] = jobs = int(self.args["--jobs"])
        if jobs == 0:
            self.args["--jobs"] = multiprocessing.cpu_count()

        return None

    def __call__(self) -> int:  # noqa: D102
        return self._kfold() if self.args["kfold"] else self._loto()

    def _load_sequences(self):
        self.logger.info("Loading the feature table")
        with open(self.args["--features"]) as in_:
            table = FeatureTable.load(in_)

        self.logger.info("Converting data to genes")
        gene_count = len(set(table.protein_id))
        genes = list(tqdm.tqdm(table.to_genes(), total=gene_count))
        del table

        self.logger.info("Sorting genes by location")
        genes.sort(key=operator.attrgetter("source.id", "start", "end"))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

        self.logger.info("Grouping genes by source sequence")
        groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
        seqs = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]

        if self.args["--shuffle"]:
            self.logger.info("Shuffling training data sequences")
            random.shuffle(seqs)

        return seqs

    def _loto(self) -> int:
        seqs = self._load_sequences()

        # --- CROSS-VALIDATION ------------------------------------------------
        k = 10
        splits = list(sklearn.model_selection.KFold(k).split(seqs))
        new_genes = []

        self.logger.info("Performing cross-validation")
        for i, (train_indices, test_indices) in enumerate(tqdm.tqdm(splits)):
            train_data = [gene for i in train_indices for gene in seqs[i]]
            test_data = [gene for i in test_indices for gene in seqs[i]]
            crf = ClusterCRF(
                self.args["--feature-type"],
                algorithm="lbfgs",
                overlap=self.args["--overlap"],
                c1=self.args["--c1"],
                c2=self.args["--c2"],
            )
            crf.fit(train_data, jobs=self.args["--jobs"], select=self.args["--select"])
            new_genes.extend(crf.predict_probabilities(test_data, jobs=self.args["--jobs"]))

            with open(self.args["--output"], "w" if i == 0 else "a") as out:
                FeatureTable.from_genes(new_genes).dump(out, header=i==0)

    def _kfold(self) -> int:
        seqs = self._load_sequences()

        # --- CROSS-VALIDATION ------------------------------------------------
        splits = list(LeaveOneGroupOut().split(seqs, groups=groups))
        new_genes = []

        self.logger.info("Performing cross-validation")
        for i, (train_indices, test_indices) in enumerate(tqdm.tqdm(splits)):
            train_data = [gene for i in train_indices for gene in seqs[i]]
            test_data = [gene for i in test_indices for gene in seqs[i]]
            crf = ClusterCRF(
                self.args["--feature-type"],
                algorithm="lbfgs",
                overlap=self.args["--overlap"],
                c1=self.args["--c1"],
                c2=self.args["--c2"],
            )
            crf.fit(train_data, jobs=self.args["--jobs"], select=self.args["--select"])
            new_genes.extend(crf.predict_probabilities(test_data, jobs=self.args["--jobs"]))

            with open(self.args["--output"], "w" if i == 0 else "a") as out:
                FeatureTable.from_genes(new_genes).dump(out, header=i==0)
