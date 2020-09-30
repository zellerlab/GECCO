"""Implementation of the ``gecco cv`` subcommand.
"""

import copy
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
from ...model import ClusterTable, FeatureTable
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
        seqs = self._load_sequences()

        if self.args["loto"]:
            splits = self._loto_splits(seqs)
        else:
            splits = self._kfold_splits(seqs)

        self.logger.info("Performing cross-validation")
        for i, (train_indices, test_indices) in enumerate(tqdm.tqdm(splits)):
            # extract train data
            train_data = [gene for i in train_indices for gene in seqs[i]]

            # extract test data and erase existing probabilities
            test_data = [copy.deepcopy(gene) for i in test_indices for gene in seqs[i]]
            for gene in test_data:
                for domain in gene.protein.domains:
                    domain.probability = None

            # fit and predict the CRF for the current fold
            crf = ClusterCRF(
                self.args["--feature-type"],
                algorithm="lbfgs",
                overlap=self.args["--overlap"],
                c1=self.args["--c1"],
                c2=self.args["--c2"],
            )
            crf.fit(train_data, jobs=self.args["--jobs"], select=self.args["--select"])
            new_genes = crf.predict_probabilities(test_data, jobs=self.args["--jobs"])

            with open(self.args["--output"], "a" if i else "w") as out:
                frame = FeatureTable.from_genes(new_genes).to_dataframe()
                frame.assign(fold=i).to_csv(out, header=i==0, sep="\t")

        return 0

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

    def _loto_splits(self, seqs):
        self.logger.info("Loading the clusters table")
        with open(self.args["--clusters"]) as in_:
            table = ClusterTable.load(in_)
            index = { protein: row.type for row in table for protein in row.proteins }

        groups = []
        for cluster in seqs:
            ty = next(index[seq.id] for seq in cluster if seq.id in index)
            groups.append(ty.unpack())

        return list(LeaveOneGroupOut().split(seqs, groups=groups))

    def _kfold_splits(self, seqs):
        k = self.args["--splits"]
        return list(sklearn.model_selection.KFold(k).split(seqs))
