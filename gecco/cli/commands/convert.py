"""Implementation of the ``gecco convert`` subcommand.
"""

import contextlib
import copy
import csv
import errno
import functools
import itertools
import glob
import os
import operator
import multiprocessing
import random
import re
import signal
import typing
import urllib.parse
from typing import Any, Dict, Union, Optional, List, TextIO, Mapping, Set

from ... import __version__
from ._base import Command, CommandExit, InvalidArgument
from .._utils import patch_showwarnings


class Convert(Command):  # noqa: D101

    summary = "convert output for compatibility with other tools"

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        return f"""
        gecco convert  - {cls.summary}

        Usage:
            gecco convert gbk -i <input> -f <format> [options]
            gecco convert clusters -i <input> -f <format> [options]

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

        ``gecco convert gbk --format=bigslice``
            Convert and alias the GenBank files in the given directory so that
            they can be loaded by BiG-SLiCE. Output files are named
            ``*.regionNNN.gbk``.

        ``gecco convert gbk --format=fna``
            Convert the GenBank files to FASTA files containing the nucleotide
            sequence of the cluster. Output files are named ``*.fna``.

        ``gecco convert gbk --format=faa``
            Convert the GenBank files to FASTA files containing the amino-acid
            sequences of all the proteins in a cluster. Output files are
            named ``*.faa``.

        ``gecco convert clusters --format=gff``
            Convert the clusters tables to GFF3 containing the position
            and metadata for all the predicted clusters. Output file is 
            named ``*.clusters.gff``.

        """

    _CLUSTERS_FORMATS: Set[str] = {"gff"}
    _GBK_FORMATS = {"bigslice", "faa", "fna"}

    def _check(self) -> None:
        try:
            super()._check()
            formats = self._GBK_FORMATS if self.args["gbk"] else self._CLUSTERS_FORMATS
            self.format = self._check_flag("--format", str, lambda x: x in formats, hint=", ".join(formats))
            self.input_dir: str = self._check_flag("--input-dir")
        except InvalidArgument:
            raise CommandExit(1)

    # ---

    def _convert_gbk_bigslice(self, ctx: contextlib.ExitStack) -> None:
        import Bio.SeqIO
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from ...model import ClusterTable

        # locate original folder
        if not os.path.exists(self.input_dir):
            self.error("Could not find input folder:", repr(self.input_dir))
            raise CommandExit(errno.ENOENT)
        # collect `*.clusters.tsv` files
        cluster_files = glob.glob(os.path.join(self.input_dir, "*.clusters.tsv"))
        unit = "table" if len(cluster_files) == 1 else "tables"
        task = self.progress.add_task("Loading clusters", total=len(cluster_files), unit=unit, precision="")
        # load the original coordinates from the `*.clusters.tsv` files
        coordinates = {}
        types = {}
        for cluster_file in self.progress.track(cluster_files, task_id=task):
            table = ClusterTable.load(cluster_file)
            for i, cluster_id in enumerate(table.cluster_id):
                coordinates[cluster_id] = (table.start[i], table.end[i])
                types[cluster_id] = table.type[i] or "Unknown"

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
                self.warn(f"GenBank file {gbk_file!r} was not obtained by GECCO")
                continue
            # mark as AntiSMASH v5
            record.annotations['structured_comment']['antiSMASH-Data'] = {
                "Version": "5.X",
                "Orig. start": coordinates[record.id][0],
                "Orig. end": coordinates[record.id][1],
            }
            # using a /subregion feature will make BiG-SLiCE think the
            # cluster is coming from MIBiG, and the predicted type will be
            # handled correctly
            subregion_feature = SeqFeature(FeatureLocation(0, len(record)), type="subregion")
            subregion_feature.qualifiers["contig_edge"] = ["False"]
            subregion_feature.qualifiers["aStool"] = ["mibig"]
            subregion_feature.qualifiers["label"] = [types[record.id]]
            record.features.append(subregion_feature)
            # rename the {id}_cluster_{N}.gbk file to {id}.region{N}.gbk
            gbk_name = os.path.basename(gbk_file)
            contig_id, cluster_n = re.search("^(.*)_cluster_(\d+).gbk", gbk_file).groups()  # type: ignore
            new_name = os.path.join(self.input_dir, "{}.region{:03}.gbk".format(contig_id, int(cluster_n)))
            self.info(f"Rewriting {gbk_file!r} to {new_name!r}")
            Bio.SeqIO.write(record, new_name, "genbank")
            done += 1

        self.success("Converted", f"{done} GenBank {unit} to BiG-SLiCE format", level=0)

    def _convert_gbk_fna(self, ctx: contextlib.ExitStack) -> None:
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
                self.warn(f"GenBank file {gbk_file!r} was not obtained by GECCO")
                continue
            # convert to nucleotide FASTA
            new_name = re.sub(r"\.gbk$", ".fna", gbk_file)
            self.info(f"Converting {gbk_file!r} to FASTA file {new_name!r}")
            Bio.SeqIO.write(record, new_name, "fasta")
            done += 1
        self.success("Converted", f"{done} GenBank {unit} to nucleotide FASTA format", level=0)

    def _convert_gbk_faa(self, ctx: contextlib.ExitStack) -> None:
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
                self.warn(f"GenBank file {gbk_file!r} was not obtained by GECCO")
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
            new_name = re.sub(r"\.gbk$", ".faa", gbk_file)
            self.info(f"Converting {gbk_file!r} proteins to FASTA file {new_name!r}")
            Bio.SeqIO.write(proteins, new_name, "fasta")
            done += 1
        self.success("Converted", f"{done} GenBank {unit} to protein FASTA format", level=0)

    def _convert_clusters_gff(self, ctx: contextlib.ExitStack) -> None:
        import Bio.SeqIO
        from ...model import ClusterTable

        # collect `*_clusters_{N}.gbk` files
        tsv_files = glob.glob(os.path.join(self.input_dir, "*.clusters.tsv"))
        unit = "file" if len(tsv_files) == 1 else "files"
        task = self.progress.add_task("Converting", total=len(tsv_files), unit=unit, precision="")
        done = 0
        # rewrite GenBank files
        for tsv_file in self.progress.track(tsv_files, task_id=task, total=len(tsv_files)):
            table = ClusterTable.load(tsv_file)
            gff_file = re.sub(r"\.tsv$", ".gff", tsv_file)
            with open(gff_file, "w") as dst:
                writer = csv.writer(dst, dialect="excel-tab")
                writer.writerow(["##gff-version 3"])
                for row in table.data.rows(named=True):
                    # load the GenBank files of the cluster to extract the GECCO version
                    gbk_file = os.path.join(self.input_dir, f"{row['cluster_id']}.gbk")
                    cluster = Bio.SeqIO.read(gbk_file, "genbank")
                    annotations = cluster.annotations['structured_comment']['GECCO-Data']
                    # make sure to have a BGC type to write down
                    bgc_types = ["Unknown"] if not row["type"] else row["type"].split(";")
                    # extract type probabilities
                    type_probas = []
                    for key, value in row.items():
                        if key.endswith("_probability"):
                            ty = key.split("_")[0].capitalize()
                            if ty == "Nrp":
                                ty = "NRP"
                            type_probas.append(f"Type{ty}={value}")
                    # write the GFF row
                    writer.writerow([
                        row["sequence_id"],
                        annotations["version"],
                        "BGC",
                        str(row["start"]),
                        str(row["end"]),
                        str(row["average_p"]),
                        ".",
                        ".",
                        ";".join([
                            f"ID={row['cluster_id']}",
                            f"Name={ '/'.join(sorted(bgc_types)) } cluster",
                            f"Type={ ','.join(sorted(bgc_types)) }",
                            f"ProbabilityAverage={row['average_p']}",
                            f"ProbabilityMax={row['max_p']}",
                            *type_probas,
                            f"Genes={row['proteins'].count(';') + 1}",
                            f"Domains={row['domains'].count(';') + 1}",
                        ])
                    ])
            done += 1
        self.success("Converted", f"{done} TSV {unit} to GFF format", level=0)

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check the CLI arguments were fine and enter context
            self._check()
            ctx.enter_context(self.progress)
            ctx.enter_context(patch_showwarnings(self._showwarnings))  # type: ignore
            # run the appropriate method
            if self.args["gbk"]:
                if self.args["--format"] == "bigslice":
                    self._convert_gbk_bigslice(ctx)
                elif self.args["--format"] == "fna":
                    self._convert_gbk_fna(ctx)
                elif self.args["--format"] == "faa":
                    self._convert_gbk_faa(ctx)
            elif self.args["clusters"]:
                if self.args["--format"] == "gff":
                    self._convert_clusters_gff(ctx)
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
