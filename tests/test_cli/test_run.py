import contextlib
import io
import os
import shutil
import tempfile
import textwrap
import unittest
from unittest import mock

import Bio.SeqIO
from Bio.Seq import Seq

import gecco.cli.commands.run
from gecco.model import ClusterTable, FeatureTable
from gecco.types import TypeClassifier
from gecco.cli import main
from gecco.cli.commands.run import Run

from ._base import TestCommand


class TestRun(TestCommand, unittest.TestCase):
    command_type = Run

    @property
    def folder(self):
        return os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_fasta_genome(self):
        sequence = os.path.join(self.folder, "data", "BGC0001866.fna")
        source = Bio.SeqIO.read(sequence, "fasta")
        with open(os.path.join(self.folder, "data", "BGC0001866.features.tsv")) as f:
            features = FeatureTable.load(f)
            genes = list(features.to_genes())

        # we mock time consuming operations (type prediction and HMM annotation)
        # with precomputed or fake results
        _find_genes = mock.Mock(return_value=genes)
        _run = mock.Mock()
        _fit_predict = lambda self, x: x
        _predict_probabilities = mock.Mock(return_value=genes)

        #_concat = mock.Mock(return_value=feats_df)
        with contextlib.ExitStack() as stack:
            # stack.enter_context(
            #     mock.patch.object(gecco.cli.commands.run.PyHMMER, "run", new=_run)
            # )
            # stack.enter_context(
            #     mock.patch.object(gecco.cli.commands.run.PyrodigalFinder, "find_genes", new=_find_genes)
            # )
            argv = ["-vv", "run", "--genome", sequence, "--output", self.tmpdir]
            with io.StringIO() as stderr:
                retcode = main(argv, stream=stderr)
                self.assertEqual(retcode, 0, stderr.getvalue())

        # make sure we have generated the files we want
        # and that we found one cluster
        output = os.listdir(self.tmpdir)
        self.assertIn("BGC0001866.features.tsv", output)
        self.assertIn("BGC0001866.clusters.tsv", output)

        with open(os.path.join(self.tmpdir, "BGC0001866.clusters.tsv")) as f:
            clusters = ClusterTable.load(f)
        self.assertEqual(len(clusters), 1)
