"""Test `gecco.orf` module.
"""

import os
import unittest
from unittest import mock

import Bio.SeqIO
from gecco.model import Strand
from gecco.orf import CDSFinder


class TestCDSFinder(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # - BGC0001737.gbk downloaded from MIBiG:
        # https://mibig.secondarymetabolites.org/repository/BGC0001737/BGC0001737.gbk
        folder = os.path.dirname(os.path.abspath(__file__))
        cls.genome = Bio.SeqIO.read(os.path.join(folder, "data", "BGC0001377.gbk"), "genbank")

    def test_sequence_coordinates(self):
        """Test emitted genes have the protein sequence matching their coordinates.
        """
        finder = CDSFinder()
        for gene in finder.find_genes([self.genome]):
            subseq = self.genome.seq[gene.start-1:gene.end]
            if gene.strand is Strand.Reverse:
                subseq = subseq.reverse_complement()
            # Don't compare first letter because it can be set as `M` in the
            # CDS feature but actually correspond to a different letter when
            # translating directly
            self.assertEqual(subseq.translate()[1:], gene.protein.seq[1:])

    def test_progress_callback(self):
        """Test that the progress callback is called for each genome.
        """
        progress = mock.MagicMock()
        finder = CDSFinder()
        genes = list(finder.find_genes([self.genome], progress=progress))
        progress.assert_called_with(self.genome, 32)
