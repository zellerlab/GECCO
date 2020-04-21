import functools
import itertools
import os
import multiprocessing
import random
import typing
from typing import List

import numpy
import pandas
import tqdm

from ._base import Command
from ...crf import ClusterCRF

if typing.TYPE_CHECKING:
    from pandas import DataFrame


class Tune(Command):

    summary = "optimize value of hyperparameters through cross-validation."
    doc = f"""
    gecco tune  - {summary}

    Usage:
        gecco tune (-h | --help)
        gecco tune kfold -i <data> [--c1 <n>]... [--c2 <n>]... [-w <col>]... [-f <col>]... [options]
        gecco tune loto  -i <data> [--c1 <n>]... [--c2 <n>]... [-w <col>]... [-f <col>]... [options]

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
        --c1 <n>                        the different values to use for C1. Can
                                        be given more than once to create a
                                        grid. [default: 0 0.15 1 2 10]
        --c2 <n>                        the different values to use for C2. Can
                                        be given more than once to create a
                                        grid. [default: 0 0.15]
        --feature-type <type>           how features should be extracted
                                        (single, overlap, or group).
                                        [default: group]
        --truncate <N>                  the maximum number of rows to use from
                                        the training set.
        --overlap <N>                   how much overlap to consider if
                                        features overlap. [default: 2]
        --splits <N>                    number of folds for cross-validation
                                        (if running `kfold`). [default: 10]
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
        self.args["--splits"] = int(self.args["--splits"])
        self.args["--c1"] = [float(c1) for c1 in self.args["--c1"]]
        self.args["--c2"] = [float(c2) for c2 in self.args["--c2"]]
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1

        # Check the `--jobs`flag
        self.args["--jobs"] = jobs = int(self.args["--jobs"])
        if jobs == 0:
            self.args["--jobs"] = multiprocessing.cpu_count()

        return None

    def __call__(self) -> int:
        # --- LOADING AND PREPROCESSING --------------------------------------
        self.logger.info("Loading the data from {!r}", self.args['--input'])
        data_tbl = pandas.read_csv(self.args["--input"], sep="\t", encoding="utf-8")
        self.logger.debug("Filtering results with e-value under {}", self.args["--e-filter"])
        data_tbl = data_tbl[data_tbl["i_Evalue"] < self.args["--e-filter"]]
        self.logger.debug("Splitting input by column {!r}", self.args["--split-col"])
        data: List['DataFrame'] = [s for _, s in data_tbl.groupby(self.args["--split-col"])]

        if self.args["--shuffle"]:
            self.logger.debug("Shuffling rows")
            random.shuffle(data)

        # --- CROSS VALIDATION -----------------------------------------------
        try:
            results = {}

            # compute results for all values
            grid = itertools.product(self.args["--c1"], self.args["--c2"])
            for c1, c2 in grid:
                # create a new CRF with C1/C2 parameters
                crf = ClusterCRF(
                    feature_columns = self.args["--feature-cols"],
                    weight_columns = self.args["--weight-cols"],
                    feature_type = self.args["--feature-type"],
                    label_column = self.args["--y-col"],
                    overlap = self.args["--overlap"],
                    algorithm = "lbfgs",
                    c1 = c1,
                    c2 = c2
                )
                # choose the right cross-validation method
                if self.args["loto"]:
                    cv_type = "loto"
                    cross_validate = crf.loto_cv
                elif self.args["kfold"]:
                    cv_type = "kfold"
                    cross_validate = functools.partial(crf.cv, k=self.args["--splits"])
                # run the cross validation
                self.logger.info("Performing cross-validation with C1={:.02}, C2={:.02}", c1, c2)
                raw = cross_validate(
                    data,
                    self.args["--strat-col"],
                    trunc=self.args["--truncate"],
                    jobs=self.args["--jobs"],
                )
                if raw:
                    results[c1, c2] = pandas.concat(raw).assign(c1=c1, c2=c2)

            return 0
        finally:
            # still attempt to write results if the tuning was
            # interrupted before completion
            if results:
                # Concatenate results
                table = pandas.concat(results.values())
                table["feature_type"] = self.args["--feature-type"]
                table["e_filter"] = self.args["--e-filter"]
                table["overlap"] = self.args["--overlap"]
                table["weight"] = ",".join(map(str, self.args["--weight-cols"]))
                table["feature"] = ",".join(self.args["--feature-cols"])
                table["truncate"] = self.args["--truncate"]
                table["input"] = os.path.basename(self.args["--input"])
                table["cv_type"] = cv_type

                # Write results
                self.logger.info("Writing output to {!r}", self.args["--output"])
                table.to_csv(self.args["--output"], sep="\t", index=False)
