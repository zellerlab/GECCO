import contextlib
import io
import os
import shutil
import tempfile
import textwrap
import unittest
from unittest import mock

import Bio.SeqIO
import pandas
from Bio.Seq import Seq

import gecco.cli.commands.run
from gecco.model import Domain, Gene, Protein, Strand
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
        sequence = os.path.join(self.folder, "BGC0001866.fna")
        source = Bio.SeqIO.read(sequence, "fasta")
        with open(os.path.join(self.folder, "BGC0001866.features.tsv")) as f:
            feats_df = pandas.read_table(f)

        genes = []
        for prot_id, df in feats_df.groupby("protein_id"):
            prot = Protein(prot_id, seq=Seq("M"))
            gene = Gene(source, min(df.start), max(df.end), Strand.Coding, prot)
            for t in df.itertuples():
                d = Domain(t.domain, t.domain_start, t.domain_end, t.hmm, t.i_Evalue, t.p_pred, {})
                gene.protein.domains.append(d)
            genes.append(gene)

        # we mock time consuming operations (type prediction and HMM annotation)
        # with precomputed or fake results
        _find_genes = mock.Mock(return_value=genes)
        _run = mock.Mock()
        _fit_predict = lambda self, x: x
        _predict_probabilities = mock.Mock(return_value=genes)

        #_concat = mock.Mock(return_value=feats_df)
        with contextlib.ExitStack() as stack:
            stack.enter_context(
                mock.patch.object(gecco.cli.commands.run.HMMER, "run", new=_run)
            )
            stack.enter_context(
                mock.patch.object(gecco.cli.commands.run.PyrodigalFinder, "find_genes", new=_find_genes)
            )
            stack.enter_context(
                mock.patch("gecco.crf.ClusterCRF.predict_probabilities", new=_predict_probabilities)
            )
            stack.enter_context(
                mock.patch("gecco.cli.commands.run.TypeClassifier.trained", new=lambda model: TypeClassifier())
            )
            stack.enter_context(
                mock.patch("gecco.cli.commands.run.TypeClassifier.predict_types", new=_fit_predict)
            )
            argv = ["-vv", "--traceback", "run", "--genome", sequence, "--output", self.tmpdir]
            main(argv, stream=io.StringIO())

        # make sure we have generated the files we want
        # and that we found one cluster
        output = os.listdir(self.tmpdir)
        self.assertIn("BGC0001866.features.tsv", output)
        self.assertIn("BGC0001866.clusters.tsv", output)
        clusters = pandas.read_table(os.path.join(self.tmpdir, "BGC0001866.clusters.tsv"))
        self.assertEqual(len(clusters), 1)
