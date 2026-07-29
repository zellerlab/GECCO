"""Test `gecco.orf` module.
"""

import os
import unittest
from unittest import mock

import Bio.SeqIO
from gecco.model import Strand
from gecco.orf import PyrodigalFinder, GFFFinder


class TestPyrodigalFinder(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # - BGC0001737.fna converted to FASTA from source GenBank file:
        # https://mibig.secondarymetabolites.org/repository/BGC0001737/BGC0001737.gbk
        # - BGC0001737.faa obtained by running Prodigal CLI in metagenome mode
        #   on BGC0001737.fna
        folder = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(folder, "data", "BGC0001737.faa")) as f:
            cls.proteins = list(Bio.SeqIO.parse(f, "fasta"))
        with open(os.path.join(folder, "data", "BGC0001737.fna")) as f:
            cls.genome = Bio.SeqIO.read(f, "fasta")
        cls.gff_file = os.path.join(folder, "data", "BGC0001737.gff")

    def test_consistency(self):
        gff_finder = GFFFinder(self.gff_file)
        pyrodigal_finder = PyrodigalFinder()
        gff_genes = list(gff_finder.find_genes([self.genome]))
        pyrodigal_genes = list(gff_finder.find_genes([self.genome]))
        self.assertEqual(len(gff_genes), len(pyrodigal_genes))
        for gff_gene, pyrodigal_gene in zip(gff_genes, pyrodigal_genes):
            self.assertEqual(gff_gene.id, pyrodigal_gene.id)
            self.assertEqual(gff_gene.start, pyrodigal_gene.start)
            self.assertEqual(gff_gene.end, pyrodigal_gene.end)
            self.assertEqual(gff_gene.strand, pyrodigal_gene.strand)
            self.assertEqual(gff_gene.protein, pyrodigal_gene.protein)
