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
            # check that start is leftmost
            self.assertLess(gene.start, gene.end)

    def test_progress_callback(self):
        """Test that the progress callback is called for each genome.
        """
        progress = mock.MagicMock()
        finder = CDSFinder()
        genes = list(finder.find_genes([self.genome], progress=progress))
        progress.assert_called_with(self.genome, 32)
