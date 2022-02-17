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

import gecco.orf
import gecco.hmmer
import gecco.cli.commands.convert
from gecco.model import Cluster, ProductType, ClusterTable, FeatureTable
from gecco.cli import main
from gecco.cli.commands.run import Run
from gecco.cli.commands.convert import Convert

from ._base import TestCommand


class TestConvert(TestCommand, unittest.TestCase):
    command_type = Convert

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
        with open(os.path.join(self.folder, "data", "BGC0001866.features.tsv")) as f:
            features = FeatureTable.load(f)
            genes = [gene.with_source(source) for gene in features.to_genes()]

        # create a cluster containing the whole sequence and write it down as
        # a GenBank file to the mock "output" folder
        cluster = Cluster("BGC0001866.1_cluster_1", genes, type=ProductType.Polyketide)
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
            retcode = main(argv, stream=stderr)
            self.assertEqual(retcode, 0, stderr.getvalue())

        # check that the output folder has a new `regionXXX.gbk`
        regions = list(map(os.path.basename, glob.glob(os.path.join(self.tmpdir, "*region*.gbk"))))
        self.assertEqual(regions, ["BGC0001866.1.region001.gbk"])
