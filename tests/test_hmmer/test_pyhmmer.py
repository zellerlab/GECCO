"""Test `gecco.hmmer` module.
"""

import copy
import os
import unittest
from unittest import mock

import Bio.SeqIO
from gecco.model import Strand, Protein, Gene
from gecco.hmmer import PyHMMER, HMM


class TestPyHMMER(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        folder = os.path.dirname(os.path.abspath(__file__))
        cls.hmm = HMM(
            id="Pfam",
            version="vX.Y",
            url="http://example.com",
            path=os.path.join(folder, "data", "minipfam.hmm"),
            size=10,
            relabel_with=r"s/(PF\d+).\d+/\1/",
        )
        with open(os.path.join(folder, "data", "proteins.faa")) as f:
            cls.records = list(Bio.SeqIO.parse(f, "fasta"))
        cls.proteins = [
            Protein(record.id, record.seq)
            for record in cls.records
        ]
        cls.genes = [
            Gene(record, 1, len(prot.seq)*3+1, Strand.Coding, prot)
            for record, prot in zip(cls.records, cls.proteins)
        ]

    def test_whitelist(self):
        # run normally, we should find 3 genes with domains
        pyhmmer = PyHMMER(self.hmm, 1)
        genes = pyhmmer.run(copy.deepcopy(self.genes))
        self.assertEqual(sum(1 for gene in genes if gene.protein.domains), 3)
        # whitelist only the first HMM, so only the first domain should be annotated
        pyhmmer = PyHMMER(self.hmm, 1, whitelist={"PF10417"})
        self.assertEqual(pyhmmer.whitelist, {"PF10417"})
        genes = pyhmmer.run(copy.deepcopy(self.genes))
        self.assertEqual(sum(1 for gene in genes if gene.protein.domains), 1)
