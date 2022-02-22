"""Test `gecco.model.ClusterTable`.
"""

import itertools
import os
import unittest
import warnings
from unittest import mock

import Bio.SeqIO
from gecco.model import ClusterTable, ProductType


class TestClusterTable(unittest.TestCase):

    # fmt: off
    @classmethod
    def setUpClass(cls):
        cls.header = list(ClusterTable.__annotations__.keys())
        cls.row = [
            "BGC1", "BGC1_cluster_1", "1", "100", "1.0", "1.0", "Unknown",
            "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", 
            "BGC0001866.1_1;BGC0001866.1_2", "PF00106;PF00107;TIGR04532"
        ]
        rows = ["\t".join(x) for x in (cls.header, cls.row)]
        cls.table = ClusterTable.load(rows)

    def test_load_unknown_type(self):
        """Check loading a BGC type named *Unknown* works.
        """
        self.assertEqual(len(self.table), 1)
        self.assertEqual(self.table.type[0], ProductType.Unknown)

    def test_dump_unknown_type(self):
        lines = self.table.dumps().splitlines()
        self.assertEqual(lines[1], "\t".join(self.row))
