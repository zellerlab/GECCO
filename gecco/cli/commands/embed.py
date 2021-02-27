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

from ._base import Command
from ._error import InvalidArgument, CommandExit
from .._utils import numpy_error_context, in_context, patch_showwarnings
from ...hmmer import HMMER


class Embed(Command):  # noqa: D101

    summary = "embed BGC annotations into non-BGC contigs for training."

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
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
            -o <out>, --output <out>      the prefix used for the output files
                                          (which will be ``<prefix>.features.tsv``
                                          and ``<prefix>.clusters.tsv``).
                                          [default: embedding]
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
        super()._check()
        try:
            self.skip = self._check_flag("--skip", int, lambda x: x > 0, "positive integer")
            self.min_size = self._check_flag("--min-size", int, lambda x: x >= 0, "positive or null integer")
            self.e_filter = self._check_flag("--e-filter", float, lambda x: 0 <= x <= 1, hint="real number between 0 and 1")
            self.jobs = self._check_flag("--jobs", int, lambda x: x >= 0, hint="positive or null integer")
            self.bgc = self.args["--bgc"]
            self.no_bgc = self.args["--no-bgc"]
            self.output = self.args["--output"]
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _read_table(self, path: str) -> "pandas.DataFrame":
        self.info("Reading", "table from", repr(path), level=2)
        return pandas.read_table(path, dtype={"domain": str})

    def _read_no_bgc(self):
        self.info("Reading", "non-BGC feature table")

        # Read the non-BGC table and assign the Y column to `0`
        _jobs = os.cpu_count() if not self.jobs else self.jobs
        with multiprocessing.pool.ThreadPool(_jobs) as pool:
            rows = pool.map(self._read_table, self.no_bgc)
            no_bgc_df = pandas.concat(rows).assign(bgc_probability="0")

        # sort and reshape
        self.info("Sorting", "non-BGC feature table", level=2)
        no_bgc_df.sort_values(by=["sequence_id", "start", "domain_start"], inplace=True)
        return [
            s
            for _, s in no_bgc_df.groupby("sequence_id", sort=True)
            if s.shape[0] > self.min_size
        ]

    def _read_bgc(self):
        self.info("Reading", "BGC feature table")

        # Read the BGC table, assign the Y column to `1`
        _jobs = os.cpu_count() if not self.jobs else self.jobs
        with multiprocessing.pool.ThreadPool(_jobs) as pool:
            rows = pool.map(read_table, self.bgc)
            bgc_df = pandas.concat(rows).assign(bgc_probability="1")
            bgc_df["BGC_id"] = bgc_df.protein_id.str.split("|").str[0]

        # sort and reshape
        self.info("Sorting", "non-BGC feature table", level=2)
        bgc_df.sort_values(by=["BGC_id", "start", "domain_start"], inplace=True)
        return [s for _, s in bgc_df.groupby("BGC_id", sort=True)]

    def _check_count(self, no_bgc_list, bgc_list):
        no_bgc_count, bgc_count = len(no_bgc_list) - self.skip, len(bgc_list)
        if no_bgc_count < bgc_count:
            self.warn("Not enough non-BGC sequences to embed BGCs:", no_bgc_count, bgc_count)

    def _embed(
        self,
        no_bgc: "pandas.DataFrame",
        bgc: "pandas.DataFrame"
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
        embed = embed[embed["i_Evalue"] < self.e_filter]
        # add additional columns based on info from BGC and non-BGC
        with numpy_error_context(divide="ignore"):
            embed = embed.assign(
                sequence_id=no_bgc["sequence_id"].apply(lambda x: x).values[0],
                BGC_id=bgc["BGC_id"].values[0],
                pseudo_pos=range(len(embed)),
                rev_i_Evalue=1 - embed["i_Evalue"],
                log_i_Evalue=-numpy.log10(embed["i_Evalue"]),
            )
        # return the embedding
        return embed

    def _make_embeddings(self, no_bgc_list, bgc_list):
        _jobs = os.cpu_count() if not self.jobs else self.jobs
        with multiprocessing.pool.ThreadPool(_jobs) as pool:
            it = zip(itertools.islice(no_bgc_list, self.skips, None), bgc_list)
            embeddings = pandas.concat(pool.starmap(embed, it))

        embeddings.sort_values(by=["sequence_id", "start", "domain_start"], inplace=True)
        return embeddings

    def _write_clusters(self, embeddings):
        self.info("Writing", "clusters table to file", repr(f"{self.prefix}.clusters.tsv"))
        with open(f"{self.prefix}.clusters.tsv", "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.write_row([
                "sequence_id", "bgc_id", "start", "end", "average_p", "max_p",
                "type", "alkaloid_probability", "polyketide_probability",
                "ripp_probability", "saccharide_probability",
                "terpene_probability", "nrp_probability",
                "other_probability", "proteins", "domains"
            ])
            positives = embeddings[embeddings.BGC == 1]
            for sequence_id, domains in positives.groupby("sequence_id"):
                ty = domains.BGC_type.values[0]
                writer.write_row([
                    sequence_id,
                    domains.BGC_id.values[0],
                    domains.start.min(),
                    domains.end.max(),
                    float(domains.bgc_probability.values[0]),
                    float(domains.bgc_probability.values[0]),
                    ty,
                    int("Alkaloid" in ty),
                    int("Polyketide" in ty),
                    int("RiPP" in ty),
                    int("Saccharide" in ty),
                    int("Terpene" in ty),
                    int("NRP" in ty),
                    int("Other" in ty),
                    ";".join(sorted(set(domains.protein_id))),
                    ";".join(sorted(set(domains.domain)))
                ])

    def _write_features(self, embeddings):
        self.info("Writing", "features table to file", repr(f"{self.prefix}.features.tsv"))
        hmm_mapping = dict(PF="Pfam", TI="Tigrfam", PT="Panther", SM="smCOGs", RF="Resfams")
        columns = [
            'sequence_id', 'protein_id', 'start', 'end', 'strand', 'domain',
            'hmm', 'i_evalue', 'domain_start', 'domain_end', 'bgc_probability'
        ]
        embeddings.to_csv(f"{self.prefix}.features.tsv", columns=columns, sep="\t", index=False)

    # ---

    @in_context
    def execute(self) -> int:  # noqa: D102
        try:
            no_bgc_list = self._read_no_bgc()
            bgc_list = self._read_bgc()
            self._check_count(no_bgc_list, bgc_list)
            embeddings = self._make_embeddings(no_bgc_list, bgc_list)
            self._write_features(embeddings)
            self._write_clusters(embeddings)
        except CommandExit as cexit:
            return cexit.code
        else:
            return 0
