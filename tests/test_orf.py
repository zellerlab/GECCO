"""Test `gecco.orf` members.
"""

import itertools
import os
import unittest
import warnings
from unittest import mock

import Bio.SeqIO
import pandas
from gecco.orf import ProdigalFinder


class TestProdigalFinder(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        fna = os.path.abspath(os.path.join(__file__, "..", "data", "SRR492066.fna"))
        cls.samples = list(itertools.islice(Bio.SeqIO.parse(fna, "fasta"), 5))

    # setting PATH to nothing will prevent the subprocess to locate PRODIGAL
    @mock.patch.dict(os.environ, {'PATH':''})
    def test_not_installed(self):
        self.assertRaises(RuntimeError, ProdigalFinder)

    @unittest.skipUnless(ProdigalFinder.has_binary(), "PRODIGAL not available")
    def test_actual_data(self):
        orf_finder = ProdigalFinder()
        proteins = orf_finder.find_proteins(self.samples)
        self.assertTrue(proteins, "PRODIGAL did not find any protein")

    @mock.patch("subprocess.run")
    def test_metagenome(self, _run):
        orf_finder = ProdigalFinder(metagenome=True)
        list(orf_finder.find_proteins(self.samples))
        self.assertTrue(any(
            _run.call_args[0][0][i:i+2] == ["-p", "meta"]
            for i in range(len(_run.call_args[0][0]))
        ))

        orf_finder = ProdigalFinder(metagenome=False)
        list(orf_finder.find_proteins(self.samples))
        self.assertFalse(any(
            _run.call_args[0][0][i:i+2] == ["-p", "meta"]
            for i in range(len(_run.call_args[0][0]))
        ))
