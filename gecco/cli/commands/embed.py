"""Implementation of the ``gecco embed`` subcommand.
"""

import contextlib
import csv
import itertools
import multiprocessing.pool
import os
import signal
import typing

from ._base import Command, InvalidArgument, CommandExit
from .._utils import numpy_error_context, in_context, patch_showwarnings

if typing.TYPE_CHECKING:
    import pandas


class Embed(Command):  # noqa: D101

    summary = "embed BGC annotations into non-BGC contigs for training."

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
        return f"""
        gecco embed - {cls.summary}

        Usage:
            gecco embed [--bgc <data>]... [--no-bgc <data>]... [options]

        Arguments:
            --bgc <data>                  the path to the annotation table
                                          containing BGC-only training instances.
            --no-bgc <data>               the path to the annotation table
                                          containing non-BGC training instances.

        Parameters:
            -M <list>, --mapping <list>   an arbitrary list of which BGC
                                          should go into which contig. Ignores
                                          ``--min-size`` and ``--skip`` when
                                          provided.
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
            self.skip = self._check_flag("--skip", int, lambda x: x >= 0, "positive or null integer")
            self.min_size = self._check_flag("--min-size", int, lambda x: x >= 0, "positive or null integer")
            self.e_filter = self._check_flag("--e-filter", float, lambda x: 0 <= x <= 1, hint="real number between 0 and 1")
            self.jobs = self._check_flag("--jobs", int, lambda x: x >= 0, hint="positive or null integer")
            self.mapping = self.args["--mapping"]
            self.bgc = self.args["--bgc"]
            self.no_bgc = self.args["--no-bgc"]
            self.output = self.args["--output"]
            if self.mapping is not None:
                self.min_size = 0
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _read_table(self, path: str) -> "pandas.DataFrame":
        import pandas

        self.info("Reading", "table from", repr(path), level=2)
        return pandas.read_table(path, dtype={"domain": str})

    def _read_no_bgc(self):
        import pandas

        self.info("Reading", "non-BGC features")

        # Read the non-BGC table and assign the Y column to `0`
        _jobs = os.cpu_count() if not self.jobs else self.jobs
        with multiprocessing.pool.ThreadPool(_jobs) as pool:
            rows = pool.map(self._read_table, self.no_bgc)
            no_bgc_df = pandas.concat(rows).assign(bgc_probability=0.0)

        # sort and reshape
        self.info("Sorting", "non-BGC features by genomic coordinates", level=2)
        no_bgc_df.sort_values(by=["sequence_id", "start", "domain_start"], inplace=True)
        return [
            s
            for _, s in no_bgc_df.groupby("sequence_id", sort=True)
            if len(s.protein_id.unique()) > self.min_size
        ]

    def _read_bgc(self):
        import pandas

        self.info("Reading", "BGC features")

        # Read the BGC table, assign the Y column to `1`
        _jobs = os.cpu_count() if not self.jobs else self.jobs
        with multiprocessing.pool.ThreadPool(_jobs) as pool:
            rows = pool.map(self._read_table, self.bgc)
            bgc_df = pandas.concat(rows).assign(bgc_probability=1.0)

        # sort and reshape
        self.info("Sorting", "BGC features by genomic coordinates", level=2)
        bgc_df.sort_values(by=["sequence_id", "start", "domain_start"], inplace=True)
        return [s for _, s in bgc_df.groupby("sequence_id", sort=True)]

    def _read_mapping(self):
        import pandas

        if self.mapping is not None:
            mapping = pandas.read_table(self.mapping)
            return { t.bgc_id:t.sequence_id for t in mapping.itertuples() }
        return None

    def _check_count(self, no_bgc_list, bgc_list):
        no_bgc_count, bgc_count = len(no_bgc_list) - self.skip, len(bgc_list)
        if no_bgc_count < bgc_count:
            self.warn("Not enough non-BGC sequences to embed BGCs:", no_bgc_count, bgc_count)

    def _embed(
        self,
        no_bgc: "pandas.DataFrame",
        bgc: "pandas.DataFrame",
    ) -> "pandas.DataFrame":
        import pandas
        import numpy

        by_prots = [s for _, s in no_bgc.groupby("protein_id", sort=False)]
        # cut the input in half to insert the bgc in the middle
        index_half = len(by_prots) // 2
        before, after = by_prots[:index_half], by_prots[index_half:]
        # find the position at which the BGC is being inserted
        insert_position = (before[-1].end.values[0] + after[0].start.values[0]) // 2
        bgc_length = bgc.end.max() - bgc.start.min()
        bgc = bgc.assign(
            start=bgc.start - bgc.start.min() + insert_position,
            end=bgc.end - bgc.start.min() + insert_position,
        )
        # shift all the 3' genes after the BGC
        after = [
            x.assign(start=x.start + bgc_length, end=x.end + bgc_length)
            for x in after
        ]
        # concat the embedding together and filter by e_value
        embed = pandas.concat(before + [bgc] + after, sort=False)
        embed = embed.reset_index(drop=True)
        embed = embed[embed["i_evalue"] < self.e_filter]
        # add additional columns based on info from BGC and non-BGC
        with numpy_error_context(numpy, divide="ignore"):
            bgc_id = bgc["sequence_id"].values[0]
            sequence_id = no_bgc["sequence_id"].apply(lambda x: x).values[0]
            embed = embed.assign(sequence_id=sequence_id, BGC_id=bgc_id)
        # return the embedding
        self.success("Finished", "embedding", repr(bgc_id), "into", repr(sequence_id), level=2)
        return embed

    def _make_embeddings(self, no_bgc_list, bgc_list, mapping):
        import pandas

        self.info("Embedding", len(bgc_list), "BGCs into", len(no_bgc_list), "contigs")
        _jobs = os.cpu_count() if not self.jobs else self.jobs

        unit = "BGC" if len(bgc_list) == 1 else "BGCs"
        task = self.progress.add_task("Embedding", unit=unit, total=len(bgc_list))

        if mapping is None:
            it = zip(itertools.islice(no_bgc_list, self.skip, None), bgc_list)
        else:
            no_bgc_index = {x.sequence_id.values[0]:x for x in no_bgc_list}
            it = [(no_bgc_index[mapping[ bgc.sequence_id.values[0] ]], bgc) for bgc in bgc_list]

        embeddings = pandas.concat([
            self._embed(*args)
            for args in self.progress.track(it, task_id=task, total=len(bgc_list))
        ])

        embeddings.sort_values(by=["sequence_id", "start", "domain_start"], inplace=True)
        return embeddings

    def _write_clusters(self, embeddings):
        self.info("Writing", "clusters table to file", repr(f"{self.output}.clusters.tsv"))
        with open(f"{self.output}.clusters.tsv", "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow([
                "sequence_id", "bgc_id", "start", "end", "average_p", "max_p",
                "type", "alkaloid_probability", "polyketide_probability",
                "ripp_probability", "saccharide_probability",
                "terpene_probability", "nrp_probability",
                "other_probability", "proteins", "domains"
            ])
            positives = embeddings[embeddings.bgc_probability == 1.0]
            for sequence_id, domains in positives.groupby("sequence_id"):
                # ty = domains.BGC_type.values[0]
                writer.writerow([
                    sequence_id,
                    domains.BGC_id.values[0],
                    domains.start.min(),
                    domains.end.max(),
                    domains.bgc_probability.values[0],
                    domains.bgc_probability.values[0],
                    "Unknown",
                    # int("Alkaloid" in ty),
                    # int("Polyketide" in ty),
                    # int("RiPP" in ty),
                    # int("Saccharide" in ty),
                    # int("Terpene" in ty),
                    # int("NRP" in ty),
                    # int("Other" in ty),
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    ";".join(sorted(set(domains.protein_id))),
                    ";".join(sorted(set(domains.domain)))
                ])

    def _write_features(self, embeddings):
        self.info("Writing", "features table to file", repr(f"{self.output}.features.tsv"))
        hmm_mapping = dict(PF="Pfam", TI="Tigrfam", PT="Panther", SM="smCOGs", RF="Resfams")
        columns = [
            'sequence_id', 'protein_id', 'start', 'end', 'strand', 'domain',
            'hmm', 'i_evalue', 'domain_start', 'domain_end', 'bgc_probability'
        ]
        embeddings[columns].to_csv(f"{self.output}.features.tsv", sep="\t", index=False)

    # ---

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check arguments and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # load inputs
            no_bgc_list = self._read_no_bgc()
            bgc_list = self._read_bgc()
            mapping = self._read_mapping()
            if mapping is None:
                self._check_count(no_bgc_list, bgc_list)
            # make embeddings
            embeddings = self._make_embeddings(no_bgc_list, bgc_list, mapping)
            # write outputs
            self._write_features(embeddings)
            self._write_clusters(embeddings)
        except CommandExit as cexit:
            return cexit.code
        except KeyboardInterrupt:
            self.error("Interrupted")
            return -signal.SIGINT
        except Exception as err:
            self.progress.stop()
            raise
        else:
            return 0
