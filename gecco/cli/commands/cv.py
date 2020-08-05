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


class Cv(Command):  # noqa: D101

    summary = "perform cross validation on a training set."
    doc = f"""
    gecco cv  - {summary}

    Usage:
        gecco cv (-h | --help)
        gecco cv kfold -i <data>  [-w <col>]... [-f <col>]... [options]

    Arguments:
        -i <data>, --input <data>       a domain annotation table with regions
                                        labeled as BGCs and non-BGCs.

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
        --truncate <N>                  the maximum number of rows to use from
                                        the training set.
        --overlap <N>                   how much overlap to consider if
                                        features overlap. [default: 2]
        --splits <N>                    number of folds for cross-validation
                                        (if running `kfold`). [default: 10]
        --select <N>                    fraction of most significant features
                                        to select from the training data.
        --shuffle                       enable shuffling of stratified rows.

    Parameters - Column Names:
        -y <col>, --y-col <col>         column with class label. [default: BGC]
        -w <col>, --weight-cols <col>   columns with local weights on features.
                                        [default: rev_i_Evalue]
        -f <col>, --feature-cols <col>  column to be used as features.
                                        [default: domain]
        -s <col>, --split-col <col>     column to be used for splitting into
                                        samples, i.e different sequences
                                        [default: sequence_id]
        -g <col>, --group-col <col>     column to be used for grouping features
                                        if `--feature-type` is *group*.
                                        [default: protein_id]
        --strat-col <col>               column to be used for stratifying the
                                        samples. [default: BGC_type]
    """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check the input exists
        input_ = self.args["--input"]
        if not os.path.exists(input_):
            self.logger.error("could not locate input file: {!r}", input_)
            return 1

        # Check the `--feature-type`
        type_ = self.args["--feature-type"]
        if type_ not in {"single", "overlap", "group"}:
            self.logger.error("Invalid value for `--feature-type`: {}", type_)
            return 1

        # Check value of numeric arguments
        if self.args["--truncate"] is not None:
            self.args["--truncate"] = int(self.args["--truncate"])
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

        # --- LOADING AND PREPROCESSING --------------------------------------
        # Load the table
        self.logger.info("Loading the data")
        with open(self.args["--input"]) as in_:
            table = FeatureTable.load(in_)

        # Converting table to genes and sort by location
        genes = sorted(table.to_genes(), key=operator.attrgetter("source.id", "start", "end"))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

        # group by sequence
        groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
        seqs = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]

        # shuffle if required
        if self.args["--shuffle"]:
            random.shuffle(seqs)

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


        with open(self.args["--input"]+".cv.tsv", "w") as out:
            FeatureTable.from_genes(new_genes).dump(out)
