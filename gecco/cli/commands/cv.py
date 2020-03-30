import functools
import os
import multiprocessing
import random
import typing
from typing import List

import pandas

from ._base import Command
from ..crf import ClusterCRF

if typing.TYPE_CHECKING:
    from pandas import DataFrame


class Cv(Command):

    summary = "perform cross validation on a training set."
    doc = f"""
    gecco cv  - {summary}

    Usage
        gecco cv (-h | --help)
        gecco train kfold -i <data>  [-w <col>]... [-f <col>]... [options]
        gecco train loto  -i <data>  [-w <col>]... [-f <col>]... [options]

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
        self.args["--min-orfs"] = int(self.args["--min-orfs"])
        self.args["--c1"] = float(self.args["--c1"])
        self.args["--c2"] = float(self.args["--c2"])
        self.args["--splits"] = int(self.args["--splits"])
        self.args["--neighbors"] = int(self.args["--neighbors"])
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1

        # Use default threshold value dependeing on postprocessing method
        if self.args["--threshold"] is None:
            if self.args["--postproc"] == "gecco":
                self.args["--threshold"] = 0.4
            elif self.args["--postproc"] == "antismash":
                self.args["--threshold"] = 0.6
        else:
            self.args["--threshold"] = float(self.args["--threshold"])

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
        self.logger.debug("Grouping columns by splitting column {!r}", self.args["--split-col"])
        data: List['DataFrame'] = [s for _, s in data_tbl.groupby(self.args["--split-col"])]

        if self.args["--shuffle"]:
            self.logger.debug("Shuffling rows")
            random.shuffle(data_tbl)

        # --- CROSS VALIDATION -----------------------------------------------
        crf = ClusterCRF(
            Y_col = self.args["--y-col"],
            feature_cols = self.args["--feature-col"],
            weight_cols = self.args["--weight-cols"],
            feature_type = self.args["--feature-type"],
            overlap = self.args["--overlap"],
            algorithm = "lbfgs",
            c1 = self.args["--c1"],
            c2 = self.args["--c2"]
        )

        if self.args["loto"]:
            cv_type = "loto"
            cross_validate = functools.partial(
                crf.loto_cv,
                type_col=self.args["--strat-col"],
                threads=self.args["--jobs"],
                trunc=self.args["--trunc"]
            )
        elif self.args["kfold"]:
            cv_type = "kfold"
            cross_validate = functools.partial(
                crf.cv,
                k=self.args["--splits"],
                strat_col=self.args["--strat-col"],
                threads=self.args["--jobs"],
                trunc=self.args["--trunc"]
            )
        else:
            raise ValueError("unknown cross-validation method")

        self.args.info("Starting cross-validation")
        results = pandas.concat(cross_validate(data_tbl))

        self.args.info("Formatting results")
        results["c1"] = self.args["--c1"]
        results["c2"] = self.args["--c2"]
        results["feature_type"] = self.args["--feature-type"]
        results["e_filter"] = self.args["--e-filter"]
        results["overlap"] = self.args["--overlap"]
        results["weight"] = ",".join(map(str, self.args["--weight-cols"]))
        results["feature"] = ",".join(self.args["--feature-cols"])
        results["truncate"] = self.args["--truncate"]
        results["input"] = os.path.basename(self.args["--input"])
        results["cv_type"] = cv_type

        # result_df = (pd .concat(results)
        #                 .assign(c1 = self.args["--c1"],
        #                     c2 = self.args["--c2"],
        #                     feature_type = self.args["--feature-type"],
        #                     e_filter = self.args["-"],
        #                     overlap = overlap,
        #                     weight = ",".join(map(str, weight_col)),
        #                     feature = ",".join(feature_col),
        #                     truncate = trunc,
        #                     in_file = ,
        #                     cv_round = "all",
        #                     cv_type = "10-fold")
        #                 .loc[ : , ["BGC", "BGC_id", "protein_id", "domain", "idx",
        #                     "p_pred", "c1", "c2", "feature_type", "e_filter", "overlap",
        #                     "weight", "truncate", "cv_type", "cv_round", "in_file"] ])

        # Write results
        self.info("Writing output to {!r}", self.args["--output"])
        results.to_table(self.args["--output"], index=False)
