"""Test `gecco.orf` module.
"""

import os
import unittest
from unittest import mock

import Bio.SeqIO
from gecco.model import Strand
from gecco.orf import PyrodigalFinder


class TestPyrodigalFinder(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # - BGC0001737.fna converted to FASTA from source GenBank file:
        # https://mibig.secondarymetabolites.org/repository/BGC0001737/BGC0001737.gbk
        # - BGC0001737.faa obtained by running Prodigal CLI in metagenome mode
        #   on BGC0001737.fna
        folder = os.path.dirname(os.path.abspath(__file__))
        cls.genome = Bio.SeqIO.read(os.path.join(folder, "BGC0001737.fna"), "fasta")
        cls.proteins = list(Bio.SeqIO.parse(os.path.join(folder, "BGC0001737.faa"), "fasta"))

    def test_order(self):
        """Test proteins are emitted in the right order from the source file.
        """
        finder = PyrodigalFinder(cpus=1)
        for expected, actual in zip(self.proteins, finder.find_genes([self.genome])):
            self.assertEqual(expected.id, actual.protein.id)

    def test_sequence_letters(self):
        """Test emitted proteins have the right sequence.
        """
        by_id = { expected.id: expected for expected in self.proteins }
        finder = PyrodigalFinder()
        for gene in finder.find_genes([self.genome]):
            self.assertEqual( gene.protein.seq, by_id[gene.id].seq )

    def test_sequence_coordinates(self):
        """Test emitted genes have the protein sequence matching their coordinates.
        """
        finder = PyrodigalFinder()
        for gene in finder.find_genes([self.genome]):
            subseq = self.genome.seq[gene.start-1:gene.end]
            if gene.strand is Strand.Reverse:
                subseq = subseq.reverse_complement()
            self.assertEqual(subseq.translate(), gene.protein.seq)

    def test_progress_callback(self):
        progress = mock.MagicMock()
        finder = PyrodigalFinder(cpus=1)
        genes = list(finder.find_genes([self.genome], progress=progress))

        progress.assert_called_with(self.genome, 10, 1)
