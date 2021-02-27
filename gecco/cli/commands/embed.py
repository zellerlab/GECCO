"""Implementation of the ``gecco embed`` subcommand.
"""

import csv
import itertools
import logging
import math
import multiprocessing.pool
import os
import pickle
import random
import typing
import warnings

import numpy
import pandas
import tqdm

from ._base import Command
from .._utils import numpy_error_context
from ...hmmer import HMMER


class Embed(Command):  # noqa: D101

    summary = "embed BGC annotations into non-BGC contigs for training."

    @classmethod
    def doc(cls, fast=False):
        return f"""
        gecco embed - {cls.summary}

        Usage:
            gecco embed (-h | --help)
            gecco embed [--bgc <data>]... [--no-bgc <data>]... [options]

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
            -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                          multithreading. Use 0 to use all of the
                                          available CPUs. [default: 0]
            --skip <N>                    skip the first N contigs while creating
                                          the embedding. [default: 0]

        """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check value of numeric arguments
        self.args["--skip"] = int(self.args["--skip"])
        self.args["--min-size"] = int(self.args["--min-size"])
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1

        # Check the `--jobs`flag
        self.args["--jobs"] = jobs = int(self.args["--jobs"])
        if jobs == 0:
            self.args["--jobs"] = multiprocessing.cpu_count()

        # Check the input exists
        for input_ in itertools.chain(self.args["--bgc"], self.args["--no-bgc"]):
            if not os.path.exists(input_):
                self.logger.error("could not locate input file: {!r}", input_)
                return 1

        return None

    def __call__(self) -> int:  # noqa: D102
        # Load input
        self.logger.info("Reading BGC and non-BGC feature tables")

        def read_table(path: str) -> "pandas.DataFrame":
            self.logger.debug("Reading table from {!r}", path)
            return pandas.read_table(path, dtype={"domain": str})

        # Read the non-BGC table, assign the Y column to `0`, sort and reshape
        with multiprocessing.pool.ThreadPool(self.args["--jobs"]) as pool:
            rows = pool.map(read_table, self.args["--no-bgc"])
            no_bgc_df = pandas.concat(rows).assign(BGC="0")
        self.logger.debug("Sorting non-BGC table")
        no_bgc_df.sort_values(by=["sequence_id", "start", "domain_start"], inplace=True)
        no_bgc_list = [
            s
            for _, s in no_bgc_df.groupby("sequence_id", sort=False)
            if s.shape[0] > self.args["--min-size"]
        ]

        # Read the BGC table, assign the Y column to `1`, and sort
        with multiprocessing.pool.ThreadPool(self.args["--jobs"]) as pool:
            rows = pool.map(read_table, self.args["--bgc"])
            bgc_df = pandas.concat(rows).assign(BGC="1")
            bgc_df["BGC_id"] = bgc_df.protein_id.str.split("|").str[0]
        self.logger.debug("Sorting BGC table")
        bgc_df.sort_values(by=["BGC_id", "start", "domain_start"], inplace=True)
        bgc_list = [s for _, s in bgc_df.groupby("BGC_id", sort=True)]

        # Check we have enough non-BGC contigs to fit the BGCs into
        no_bgc_count, bgc_count = len(no_bgc_list) - self.args["--skip"], len(bgc_list)
        if no_bgc_count < bgc_count:
            msg = "Not enough non-BGC sequences to fit the BGCS: {} / {}"
            warnings.warn(msg.format(no_bgc_count, bgc_count))

        # Make a progress bar if we are printing to a terminal
        self.logger.info("Creating the embeddings")
        if self.stream.isatty() and self.logger.level != 0:
            pbar = tqdm.tqdm(total=min(len(no_bgc_list), len(bgc_list)), leave=False)
        else:
            pbar = None

        # Make the embeddings
        def embed(
            no_bgc: "pandas.DataFrame", bgc: "pandas.DataFrame"
        ) -> "pandas.DataFrame":
            by_prots = [s for _, s in no_bgc.groupby("protein_id", sort=False)]
            # cut the input in half to insert the bgc in the middle
            index_half = len(by_prots) // 2
            before, after = by_prots[:index_half], by_prots[index_half:]
            # compute offsets
            insert_position = (before[-1].end.max() + after[0].start.min()) // 2
            bgc_length = bgc.end.max() - bgc.start.min()
            # update offsets
            bgc = bgc.assign(
                start=bgc.start + insert_position, end=bgc.end + insert_position
            )
            after = [
                x.assign(start=x.start + bgc_length, end=x.end + bgc_length)
                for x in after
            ]
            # concat the embedding together and filter by e_value
            embed = pandas.concat(before + [bgc] + after, sort=False)
            embed = embed.reset_index(drop=True)
            embed = embed[embed["i_Evalue"] < self.args["--e-filter"]]
            # add additional columns based on info from BGC and non-BGC
            with numpy_error_context(divide="ignore"):
                embed = embed.assign(
                    sequence_id=no_bgc["sequence_id"].apply(lambda x: x).values[0],
                    BGC_id=bgc["BGC_id"].values[0],
                    pseudo_pos=range(len(embed)),
                    rev_i_Evalue=1 - embed["i_Evalue"],
                    log_i_Evalue=-numpy.log10(embed["i_Evalue"]),
                )
            # Update the progressbar, if any
            if pbar is not None:
                pbar.update(1)
            return embed

        with multiprocessing.pool.ThreadPool(self.args["--jobs"]) as pool:
            it = zip(itertools.islice(no_bgc_list, self.args["--skip"], None), bgc_list)
            embeddings = pandas.concat(pool.starmap(embed, it))
        if pbar is not None:
            pbar.close()

        # Write the resulting table
        embeddings.sort_values(
            by=["sequence_id", "start", "domain_start"], inplace=True
        )
        self.logger.info("Writing embedding table")
        out_file = self.args["--output"]
        self.logger.debug("Writing embedding table to {!r}", out_file)
        embeddings.to_csv(out_file, sep="\t", index=False)
        return 0
