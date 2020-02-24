import csv
import logging
import multiprocessing
import os
import pickle
import random
import typing

import numpy
import pandas
from Bio import SeqIO

from ._base import Command
from ... import data
from ...crf import ClusterCRF
from ...hmmer import HMMER
from ...knn import ClusterKNN
from ...orf import ORFFinder
from ...refine import ClusterRefiner
from ...preprocessing import truncate


class Train(Command):

    summary = "train the CRF model on an embedded feature table."
    doc = f"""
    gecco train - {summary}

    Usage:
        gecco train (-h | --help)
        gecco train -i <data> [-w <col>]... [--sort-cols <col>]...
                    [--strat-cols <col>]... [options]

    Arguments:
        -i <data>, --input <data>      a FASTA or GenBank file containing a
                                       genome as input.
    Parameters:
        -o <out>, --output <out>       the basename to use for the output
                                       model. [default: CRF]
        -j <jobs>, --jobs <jobs>       the number of CPUs to use for
                                       multithreading. Use 0 to use all of the
                                       available CPUs. [default: 0]
        -e <e>, --e-filter <e>         the e-value cutoff for domains to
                                       be included [default: 1e-5]
        -m <m>, --threshold <m>        the probability threshold for cluster
                                       prediction. Default depends on the
                                       post-processing method (0.4 for gecco,
                                       0.6 for antismash).
        --c1 <C1>                      parameter for L1 regularisation.
                                       [default: 0.15]
        --c2 <C2>                      parameter for L2 regularisation.
                                       [default: 0.15]
        --feature-type <type>          how features should be extracted
                                       (single, overlap, or group).
                                       [default: group]
        --truncate <N>                 the maximum number of rows to use from
                                       the training set.
        --overlap <N>                  how much overlap to consider if
                                       features overlap. [default: 2]
        --min-orfs <N>                 how many ORFs are required for a
                                       sequence to be considered. [default: 5]
        --postproc <method>            the method used for cluster extraction
                                       (antismash or gecco). [default: gecco]

    Parameters - Column Names:
        -y <col>, --y-col <col>        column with class label. [default: BGC]
        -w <col>, --weight-cols <col>  columns with local weights on features.
                                       [default: 1]
        -f <col>, --feature-col <col>  column to be used as features.
                                       [default: domain]
        -s <col>, --split-col <col>    column to be used for splitting into
                                       samples, i.e different sequences
                                       [default: sequence_id]
        -g <col>, --group-col <col>    column to be used for grouping features
                                       if `--feature-type` is *group*.
                                       [default: protein_id]
        --sort-cols <col>              columns to be used for sorting the data
                                       [default: genome_id start domain_start]
        --strat-cols <col>             columns to be used for stratifying the
                                       samples (BGC types).

    Parameters - Cross-Validation:
        --no-shuffle                   disable shuffling of the data before
                                       doing the cross-validation.
        --folds <N>                    the number of folds to use for the
                                       cross-validation. [default: 10]

    Parameters - Type Prediction:
        -d <d>, --distance <d>         the distance metric to use for kNN type
                                       prediction. [default: jensenshannon]
        -k <n>, --neighbors <n>        the number of neighbors to use for
                                       kNN type prediction [default: 5]
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
        self.args["--folds"] = int(self.args["--folds"])
        self.args["--overlap"] = int(self.args["--overlap"])
        self.args["--min-orfs"] = int(self.args["--min-orfs"])
        self.args["--c1"] = float(self.args["--c1"])
        self.args["--c2"] = float(self.args["--c2"])
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

        # Retype the --weight-cols
        self.args["--weight-cols"] = [int(w) for w in self.args["--weight-cols"]]

        # Check the `--jobs`flag
        self.args["--jobs"] = jobs = int(self.args["--jobs"])
        if jobs == 0:
            self.args["--jobs"] = multiprocessing.cpu_count()

        return None

    def __call__(self) -> int:
        # --- LOADING AND PREPROCESSING --------------------------------------
        # Load the table
        self.logger.info("Loading the data")
        data_tbl = pandas.read_csv(self.args["--input"], sep="\t", encoding="utf-8")
        self.logger.debug("Filtering results with e-value under {}", self.args["--e-filter"])
        data_tbl = data_tbl[data_tbl["i_Evalue"] < self.args["--e-filter"]]
        self.logger.debug("Reformating PFAM ids in column {}", self.args["--feature-col"])
        data_tbl = data_tbl.assign(
            domain = data_tbl[self.args["--feature-col"]]
                .str
                .replace(r"(PF\d+)\.\d+", lambda m: m.group(1))
        )

        # Computing reverse i_Evalue
        self.logger.debug("Computing reverse i_Evalue")
        data_tbl = data_tbl.assign(rev_i_Evalue = 1 - data_tbl["i_Evalue"])

        # Grouping column
        self.logger.debug("Splitting data using column {}", self.args["--split-col"])
        data_tbl = [s for _, s in data_tbl.groupby(self.args["--split-col"])]
        if not self.args["--no-shuffle"]:
            self.logger.debug("Shuffling data")
            random.shuffle(data_tbl)

        # Truncate the input data if required
        if self.args["--truncate"] is not None:
            self.logger.debug("Truncating data to {} rows", self.args["--truncate"])
            data_tbl = [
                truncate(
                    df,
                    self.args["--truncate"],
                    Y_col=self.args["--y-col"],
                    grouping=self.args["--group-col"]
                )
                for df in data_tbl
            ]

        # --- MODEL FITTING --------------------------------------------------
        crf = ClusterCRF(
            Y_col = self.args["--y-col"],
            feature_cols = self.args["--feature-col"],
            group_col = self.args["--group-col"],
            #weight_cols = self.args["--weight-cols"],
            weight_cols = ["rev_i_Evalue"],
            feature_type = self.args["--feature-type"],
            overlap = self.args["--overlap"],
            algorithm = "lbfgs",
            c1 = self.args["--c1"],
            c2 = self.args["--c2"]
        )
        self.logger.info("Fitting the model")
        crf.fit(data=data_tbl)

        model_out = f"{self.args['--output']}.crf.model"
        self.logger.info("Writing the model to {!r}", model_out)
        with open(model_out, "wb") as f:
            pickle.dump(crf, f, protocol=3)

        self.logger.info("Writing weights to {0}.trans.tsv and {0}.state.tsv", self.args["--output"])
        crf.save_weights(self.args["--output"])

        return 0
