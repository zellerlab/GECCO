import contextlib
import glob
import io
import os
import shutil
import tempfile
import textwrap
import unittest
from unittest import mock

import Bio.SeqIO
from Bio.Seq import Seq
from rich.console import Console

import gecco.orf
import gecco.hmmer
import gecco.cli.commands.convert
from gecco.model import Cluster, ClusterType, ClusterTable, FeatureTable
from gecco.cli import main

from ._base import TestCommand


class TestConvert(TestCommand, unittest.TestCase):
    name = "convert"

    @property
    def folder(self):
        return os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_convert_gbk_bigslice(self):
        # load the input sequence and pre-compute feature table
        sequence = os.path.join(self.folder, "data", "BGC0001866.fna")
        source = Bio.SeqIO.read(sequence, "fasta")
        path = os.path.join(self.folder, "data", "BGC0001866.features.tsv")
        features = FeatureTable.load(path)
        genes = [gene.with_source(source) for gene in features.to_genes()]

        # create a cluster containing the whole sequence and write it down as
        # a GenBank file to the mock "output" folder
        cluster = Cluster("BGC0001866.1_cluster_1", genes, type=ClusterType("Polyketide"))
        record = cluster.to_seq_record()
        with open(os.path.join(self.tmpdir, "{}.gbk".format(cluster.id)), "w") as f:
            Bio.SeqIO.write(record, f, "genbank")
        # copy the `.clusters.tsv` and `.features.tsv` files
        for tsv in ["clusters.tsv", "features.tsv"]:
            shutil.copy(
                os.path.join(self.folder, "data", "BGC0001866.{}".format(tsv)),
                os.path.join(self.tmpdir, "BGC0001866.1.{}".format(tsv))
            )

        # run the `gecco convert gbk` command
        argv = ["-vv", "convert", "gbk", "--input-dir", self.tmpdir, "--format", "bigslice"]
        with io.StringIO() as stderr:
            retcode = main(argv, console=Console(file=stderr))
            self.assertEqual(retcode, 0, stderr.getvalue())

        # check that the output folder has a new `regionXXX.gbk`
        regions = list(map(os.path.basename, glob.glob(os.path.join(self.tmpdir, "*region*.gbk"))))
        self.assertEqual(regions, ["BGC0001866.1.region001.gbk"])
