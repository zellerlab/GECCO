"""Implementation of the ``gecco convert`` subcommand.
"""

import contextlib
import copy
import errno
import functools
import itertools
import glob
import os
import operator
import multiprocessing
import random
import typing
from typing import Any, Dict, Union, Optional, List, TextIO, Mapping

from ... import __version__
from ...model import ClusterTable
from ._base import Command, CommandExit, InvalidArgument
from .._utils import patch_showwarnings


class Convert(Command):  # noqa: D101

    summary = "convert output for compatibility with other tools"

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
        return f"""
        gecco convert  - {cls.summary}

        Usage:
            gecco convert gbk -i <input> -f <format> [options]

        Arguments:
            -i <input>, --input-dir <input> the path to the input directory
                                            containing files to convert.
            -f <format>, --format <format>  the name of the output format to
                                            write.

        Parameters:
            -o <out>, --output <out>        the name of the output file,
                                            default depends on the output
                                            format.

        This command helps transforming the output files created by
        GECCO into helpful format, should you want to use the results in
        combination with other tools. The supported formats are listed
        below:

        * ``gecco convert gbk --format=bigslice``: convert and alias the
          GenBank files in the given directory so that they can be loaded by
          BiG-SLiCE. Output files are named ``*.regionNNN.gbk``.
        * ``gecco convert gbk --format=fna``: convert the GenBank files to
          FASTA files containing the nucleotide sequence of the cluster.
          Output files are named ``*.fna``.
        * ``gecco convert gbk --format=faa``: convert the GenBank files to
          FASTA files containing the amino-acid sequences of all the proteins
          in a cluster. Output files are named ``*.faa``.

        """

    _CLUSTERS_FORMATS = set()
    _GBK_FORMATS = {"bigslice", "faa", "fna"}

    def _check(self) -> typing.Optional[int]:
        try:
            super()._check()
            formats = self._GBK_FORMATS if self.args["gbk"] else self._CLUSTERS_FORMATS
            self.format = self._check_flag("--format", str, lambda x: x in formats, hint=", ".join(formats))
            self.input_dir = self._check_flag("--input-dir")
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _convert_gbk_bigslice(self, ctx: contextlib.ExitStack):
        import Bio.SeqIO
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        # locate original folder
        if not os.path.exists(self.input_dir):
            self.error("Could not find input folder:", repr(self.input_dir))
            raise CommandExit(errno.ENOENT)
        # collect `*.clusters.tsv` files
        cluster_files = glob.glob(os.path.join(self.input_dir, "*.clusters.tsv"))
        unit = "table" if len(cluster_files) == 1 else "tables"
        task = self.progress.add_task("Loading", total=len(cluster_files), unit=unit, precision="")
        # load the original coordinates from the `*.clusters.tsv` files
        coordinates = {}
        types = {}
        for cluster_file in self.progress.track(cluster_files, task_id=task, precision=""):
            cluster_fh = ctx.enter_context(open(cluster_file))
            for row in ClusterTable.load(cluster_fh):
                ty = ";".join(sorted(ty.name for ty in row.type.unpack()))
                coordinates[row.bgc_id] = (row.start, row.end)
                types[row.bgc_id] = ty or "Unknown"

        # collect `*_clusters_{N}.gbk` files
        gbk_files = glob.glob(os.path.join(self.input_dir, "*_cluster_*.gbk"))
        unit = "file" if len(gbk_files) == 1 else "files"
        task = self.progress.add_task("Converting", total=len(gbk_files), unit=unit, precision="")
        done = 0
        # rewrite GenBank files
        for gbk_file in self.progress.track(gbk_files, task_id=task, total=len(gbk_files), precision=""):
            # load record and ensure it comes from GECCO
            record = Bio.SeqIO.read(gbk_file, "genbank")
            if "GECCO-Data" not in record.annotations.get('structured_comment', {}):
                self.warning(f"GenBank file {gbk_file!r} was not obtained by GECCO")
                continue
            # mark as AntiSMASH v5
            record.annotations['structured_comment']['antiSMASH-Data'] = {
                "Version": "5.X",
                "Orig. start": coordinates[record.id][0],
                "Orig. end": coordinates[record.id][1],
            }
            # using a /subregion feature will make BiG-SLiCE think the
            # BGC is coming from MIBiG, and the predicted type will be
            # handled correctly
            subregion_feature = SeqFeature(FeatureLocation(0, len(record)), type="subregion")
            subregion_feature.qualifiers["contig_edge"] = ["False"]
            subregion_feature.qualifiers["aStool"] = ["mibig"]
            subregion_feature.qualifiers["label"] = [types[record.id]]
            record.features.append(subregion_feature)
            # rename the {id}_cluster_{N}.gbk file to {id}.region{N}.gbk
            new_name = gbk_file.replace("_cluster_", ".region")
            self.info(f"Rewriting {gbk_file!r} to {new_name!r}")
            Bio.SeqIO.write(record, new_name, "genbank")
            done += 1

        self.success(f"Converted {done} GenBank {unit} to BiG-SLiCE format", level=0)

    def _convert_gbk_fna(self, ctx: contextlib.ExitStack):
        import Bio.SeqIO
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        # collect `*_clusters_{N}.gbk` files
        gbk_files = glob.glob(os.path.join(self.input_dir, "*_cluster_*.gbk"))
        unit = "file" if len(gbk_files) == 1 else "files"
        task = self.progress.add_task("Converting", total=len(gbk_files), unit=unit, precision="")
        done = 0
        # rewrite GenBank files
        for gbk_file in self.progress.track(gbk_files, task_id=task, total=len(gbk_files)):
            # load record and ensure it comes from GECCO
            record = Bio.SeqIO.read(gbk_file, "genbank")
            if "GECCO-Data" not in record.annotations.get('structured_comment', {}):
                self.warning(f"GenBank file {gbk_file!r} was not obtained by GECCO")
                continue
            # convert to nucleotide FASTA
            new_name = gbk_file.replace(".gbk", ".fna")
            self.info(f"Converting {gbk_file!r} to FASTA file {new_name!r}")
            Bio.SeqIO.write(record, new_name, "fasta")
            done += 1
        self.success(f"Converted {done} GenBank {unit} to nucleotide FASTA format", level=0)

    def _convert_gbk_faa(self, ctx: contextlib.ExitStack):
        import Bio.SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        # collect `*_clusters_{N}.gbk` files
        gbk_files = glob.glob(os.path.join(self.input_dir, "*_cluster_*.gbk"))
        unit = "file" if len(gbk_files) == 1 else "files"
        task = self.progress.add_task("Converting", total=len(gbk_files), unit=unit, precision="")
        done = 0
        # rewrite GenBank files
        for gbk_file in self.progress.track(gbk_files, task_id=task, total=len(gbk_files)):
            # load record and ensure it comes from GECCO
            record = Bio.SeqIO.read(gbk_file, "genbank")
            if "GECCO-Data" not in record.annotations.get('structured_comment', {}):
                self.warning(f"GenBank file {gbk_file!r} was not obtained by GECCO")
                continue
            # extract proteins
            proteins = []
            for feat in filter(lambda f: f.type == "CDS", record.features):
                if "locus_tag" not in feat.qualifiers:
                    continue
                seq = Seq(feat.qualifiers["translation"][0])
                prot = SeqRecord(id=feat.qualifiers["locus_tag"][0], seq=seq)
                proteins.append(prot)
            # write proteins to FASTA
            new_name = gbk_file.replace(".gbk", ".faa")
            self.info(f"Converting {gbk_file!r} proteins to FASTA file {new_name!r}")
            Bio.SeqIO.write(proteins, new_name, "fasta")
            done += 1
        self.success(f"Converted {done} GenBank {unit} to protein FASTA format", level=0)

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check the CLI arguments were fine and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # run the appropriate method
            if self.args["gbk"]:
                if self.args["--format"] == "bigslice":
                    self._convert_gbk_bigslice(ctx)
                elif self.args["--format"] == "fna":
                    self._convert_gbk_fna(ctx)
                elif self.args["--format"] == "faa":
                    self._convert_gbk_faa(ctx)
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
