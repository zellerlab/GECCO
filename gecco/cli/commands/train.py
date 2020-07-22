"""Implementation of the ``gecco train`` subcommand.
"""

import csv
import hashlib
import itertools
import logging
import multiprocessing.pool
import os
import pickle
import random
import typing

import numpy
import pandas
import scipy.sparse
import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ._base import Command
from ...model import Domain, Gene, Protein, Strand
from ...crf import ClusterCRF
from ...refine import ClusterRefiner
from ...hmmer import HMMER


class Train(Command):  # noqa: D101

    summary = "train the CRF model on an embedded feature table."
    doc = f"""
    gecco train - {summary}

    Usage:
        gecco train (-h | --help)
        gecco train -i <data> [-w <col>]... [--feature-cols <col>]...
                    [--sort-cols <col>]... [--strat-cols <col>]... [options]

    Arguments:
        -i <data>, --input <data>       a domain annotation table with regions
                                        labeled as BGCs and non-BGCs.

    Parameters:
        -o <out>, --output-dir <out>    the directory to use for the model
                                        files. [default: CRF]
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
        --no-shuffle                    disable shuffling of the data before
                                        fitting the model.
        --select <N>                    fraction of most significant features
                                        to select from the training data.

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
        --sort-cols <col>               columns to be used for sorting the data
                                        [default: genome_id start domain_start]
        --strat-cols <col>              columns to be used for stratifying the
                                        samples (BGC types).

    Parameters - Type Prediction:
        --type-col <col>                column containing BGC types to use for
                                        domain composition. [default: BGC_type]
        --id-col <col>                  column containing BGC id to use for
                                        BGC labelling. [default: BGC_id]
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
        feats_df = pandas.read_csv(self.args["--input"], sep="\t", encoding="utf-8")
        self.logger.debug(
            "Filtering results with e-value under {}", self.args["--e-filter"]
        )
        feats_df = feats_df[feats_df["i_Evalue"] < self.args["--e-filter"]]

        # Computing reverse i_Evalue
        self.logger.debug("Computing reverse indepenent e-value")
        feats_df = feats_df.assign(rev_i_Evalue=1 - feats_df["i_Evalue"])

        # Sorting data
        self.logger.debug("Sorting data by gene and domain coordinates")
        feats_df.sort_values(
            by=[self.args["--split-col"], "start", "end", "domain_start"], inplace=True
        )

        # Grouping column
        self.logger.debug(
            "Splitting feature table using column {!r}", self.args["--split-col"]
        )
        training_data = [s for _, s in feats_df.groupby(self.args["--split-col"])]
        if not self.args["--no-shuffle"]:
            self.logger.debug("Shuffling splits randomly")
            random.shuffle(training_data)

        # --- MODEL FITTING --------------------------------------------------
        crf = ClusterCRF(
            label_column=self.args["--y-col"],
            feature_columns=self.args["--feature-cols"],
            weight_columns=self.args["--weight-cols"],
            group_column=self.args["--group-col"],
            feature_type=self.args["--feature-type"],
            overlap=self.args["--overlap"],
            algorithm="lbfgs",
            c1=self.args["--c1"],
            c2=self.args["--c2"],
        )
        self.logger.info("Fitting the CRF model to the training data")
        crf.fit(
            data=training_data,
            trunc=self.args["--truncate"],
            select=self.args["--select"],
        )

        os.makedirs(self.args["--output-dir"], exist_ok=True)
        model_out = os.path.join(self.args["--output-dir"], "model.pkl")
        self.logger.info("Writing the model to {!r}", model_out)
        with open(model_out, "wb") as f:
            pickle.dump(crf, f, protocol=4)

        self.logger.debug("Computing and saving model checksum")
        hasher = hashlib.md5()
        with open(model_out, "rb") as f:
            hasher.update(f.read())  # FIXME: iterate on file blocks
        with open(f"{model_out}.md5", "w") as f:
            f.write(hasher.hexdigest())

        self.logger.info("Writing transitions and state weights")
        crf.save_weights(self.args["--output-dir"])

        # # --- DOMAIN COMPOSITION ----------------------------------------------
        self.logger.info("Extracting domain composition for type classifier")

        self.logger.debug("Finding the array of possible domains")
        if crf.significant_features:
            all_possible = sorted(
                {d for domains in crf.significant_features.values() for d in domains}
            )
        else:
            all_possible = sorted(
                {t.domain for d in training_data for t in d[d["BGC"] == 1].itertuples()}
            )
        types = [d[self.args["--type-col"]].values[0] for d in training_data]
        ids = [d[self.args["--id-col"]].values[0] for d in training_data]

        self.logger.debug("Saving training matrix labels for BGC type classifier")
        doms_out = os.path.join(self.args["--output-dir"], "domains.tsv")
        pandas.Series(all_possible).to_csv(
            doms_out, sep="\t", index=False, header=False
        )
        types_out = os.path.join(self.args["--output-dir"], "types.tsv")
        df = pandas.DataFrame(dict(ids=ids, types=types))
        df.to_csv(types_out, sep="\t", index=False, header=False)

        self.logger.debug("Building new domain composition matrix")
        if self.stream.isatty() and self.logger.level != 0:
            pbar = tqdm.tqdm(total=len(training_data), leave=False)
        else:
            pbar = None

        def domain_composition(table: pandas.DataFrame) -> numpy.ndarray:
            is_bgc = table[self.args["--y-col"]].array == 1
            names = table[self.args["--feature-cols"][0]].array[is_bgc]
            weights = table[self.args["--weight-cols"][0]].array[is_bgc]
            composition = numpy.zeros(len(all_possible))
            for i, domain in enumerate(all_possible):
                composition[i] = numpy.sum(weights[names == domain])
            if pbar is not None:
                pbar.update(1)
            return composition / (composition.sum() or 1)  # type: ignore

        with multiprocessing.pool.ThreadPool(self.args["--jobs"]) as pool:
            new_comp = numpy.array(pool.map(domain_composition, training_data))
        if pbar is not None:
            pbar.close()

        comp_out = os.path.join(self.args["--output-dir"], "compositions.npz")
        self.logger.debug("Saving new domain composition matrix to {!r}", comp_out)
        scipy.sparse.save_npz(comp_out, scipy.sparse.coo_matrix(new_comp))
        return 0
