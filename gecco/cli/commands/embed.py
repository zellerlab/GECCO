import csv
import logging
import math
import multiprocessing
import os
import pickle
import random
import typing
import warnings

import numpy
import pandas
import tqdm
from Bio import SeqIO

from ._base import Command
from .._utils import numpy_error_context
from ... import data
from ...crf import ClusterCRF
from ...hmmer import HMMER
from ...knn import ClusterKNN
from ...orf import ORFFinder
from ...refine import ClusterRefiner
from ...preprocessing import truncate


class Embed(Command):

    summary = "embed BGC annotations into non-BGC contigs for training."
    doc = f"""
    gecco embed - {summary}

    Usage:
        gecco embed (-h | --help)
        gecco embed --bgc <data> --no-bgc <data> [options]

    Arguments:
        --bgc <data>                  the path to the annotation table
                                      containing BGC-only training instances.
        --no-bgc <data>               the path to the annotation table
                                      containing non-BGC training instances.

    Parameters:
        -o <out>, --output <out>      the file in which to write the
                                      resulting embedding table.
                                      [default: features.tsv]
        --min-size <N>                the minimum size for padding sequences.
                                      [default: 500]
        -e <e>, --e-filter <e>        the e-value cutoff for domains to be
                                      included. [default: 1e-5]

    """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check value of numeric arguments
        self.args["--min-size"] = int(self.args["--min-size"])
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1

        # Check the input exists
        for input_ in (self.args["--bgc"], self.args["--no-bgc"]):
            if not os.path.exists(input_):
                self.logger.error("could not locate input file: {!r}", input_)
                return 1

        return None

    def __call__(self) -> int:
        # Check CLI arguments
        retcode = self._check()
        if retcode is not None:
            return retcode

        self.logger.info("Reading BGC and non-BGC feature tables")

        # Read the non-BGC table, assign the Y column to `0`, sort and reshape
        self.logger.debug("Reading non-BGC table from {!r}", self.args["--no-bgc"])
        no_bgc_df = pandas.read_csv(self.args["--no-bgc"], sep="\t")
        no_bgc_df = no_bgc_df.assign(BGC="0")
        self.logger.debug("Sorting non-BGC table")
        no_bgc_df = no_bgc_df.sort_values(by=["sequence_id", "start", "domain_start"])
        no_bgc_df = no_bgc_df.groupby("sequence_id", sort=False)
        no_bgc_list = [s for _, s in no_bgc_df if s.shape[0] > self.args["--min-size"]]

        # Read the BGC table, assign the Y column to `0`, and sort
        self.logger.debug("Reading BGC table from {!r}", self.args["--bgc"])
        bgc_df = pandas.read_csv(self.args["--bgc"], sep="\t")
        bgc_df = bgc_df.assign(
            BGC="1",
            BGC_id=[id_[0] for id_ in bgc_df['protein_id'].str.split("|")]
        )
        self.logger.debug("Sorting BGC table")
        bgc_df = bgc_df.sort_values(by=["BGC_id", "start", "domain_start"])
        bgc_list = [s for _, s in bgc_df.groupby("BGC_id", sort=True)]

        # Checking we have enough non-BGC contigs to fit the BGCs into
        no_bgc_count, bgc_count = len(no_bgc_count), len(bgc_count)
        if no_bgc_count < bgc_count:
            msg = "Not enough non-BGC sequences to fit the BGCS: {} / {}"
            warnings.warn(msg.format(no_bgc_count, bgc_count))
            embedding_count = bgc_count
        else:
            embedding_count = no_bgc_count

        # Make the embeddings
        self.logger.info("Creating the embeddings")
        embedding = []
        for no_bgc, bgc in tqdm.tqdm(zip(no_bgc_list, bgc_list), total=embedding_count):
            no_bgc = no_bgc.groupby("protein_id", sort=False)
            start, end = pandas.DataFrame(), pandas.DataFrame()
            for n, (_, t) in enumerate(tqdm.tqdm(no_bgc, leave=False)):
                if n <= math.ceil(len(no_bgc) / 2):
                    start = pandas.concat([start, t])
                else:
                    end = pandas.concat([end, t])

            embed = pandas.concat([start, bgc, end], sort=False)
            embed = embed.reset_index(drop=True)
            embed = embed[embed["i_Evalue"] < self.args["--e-filter"]]

            #
            with numpy_error_context(divide="ignore"):
                embed = embed.assign(
                    domain=embed["domain"].str.replace(r"(PF\d+)\.\d+", lambda m: m.group(1)),
                    sequence_id=no_bgc["sequence_id"].apply(lambda x: x).values[0],
                    BGC_id=bgc["BGC_id"].values[0],
                    # FIXME: really needed ? if so we must also extract the
                    # metadata from the `mibig.json` metadata listing
                    # BGC_type=bgc["BGC_type"].values[0].split(","),
                    pseudo_pos=range(len(embed)),
                    rev_i_Evalue = 1 - embed["i_Evalue"],
                    log_i_Evalue = -numpy.log10(embed["i_Evalue"]),
                    strand_shift = ~embed["strand"].eq(embed["strand"].shift(1)),
                    shift = "shift"
                )
                embedding.append(embed)

        # Write the resulting table
        self.logger.info("Writing embedding table")
        out_file = self.args["--output"]
        self.logger.debug("Writing embedding table to {!r}", out_file)
        pandas.concat(embedding).to_csv(out_file, sep="\t", index=False)
