"""Test `gecco.model.ClusterTable`.
"""

import itertools
import os
import unittest
import warnings
from unittest import mock

import Bio.SeqIO
from gecco.model import ClusterTable, ClusterType


class TestClusterTable(unittest.TestCase):

    # fmt: off
    @classmethod
    def setUpClass(cls):
        cls.header = [
            "sequence_id", "cluster_id", "start", "end", "average_p",
            "max_p", "type", "alkaloid_probability", "nrp_probability",
            "polyketide_probability", "ripp_probability", "saccharide_probability",
            "terpene_probability", "proteins", "domains",
        ]
        cls.row = [
            "BGC1", "BGC1_cluster_1", "1", "100", "1.0", "1.0", "Polyketide",
            "0.6", "0.5", "0.4", "0.3", "0.2", "0.1",
            "BGC0001866.1_1;BGC0001866.1_2", "PF00106;PF00107;TIGR04532"
        ]
        rows = ["\t".join(x) for x in (cls.header, cls.row)]
        data = "\n".join(rows).encode("utf-8")
        cls.table = ClusterTable.loads(data)

    def test_load(self):
        self.assertEqual(len(self.table), 1)
        self.assertEqual(self.table.type[0], "Polyketide")

    def test_load_unknown_product_type(self):
        row = self.row.copy()
        row[6] = "Unknown"
        rows = ["\t".join(x) for x in (self.header, row)]
        data = "\n".join(rows).encode("utf-8")
        table = ClusterTable.loads(data)
        self.assertEqual(len(table), 1)
        self.assertEqual(table.type[0], "Unknown")

    def test_load_no_product_type(self):
        row = self.row.copy()
        row[6] = ""
        rows = ["\t".join(x) for x in (self.header, row)]
        data = "\n".join(rows).encode("utf-8")
        table = ClusterTable.loads(data)
        self.assertEqual(len(table), 1)
        self.assertIs(table.type[0], None)

    def test_dump(self):
        lines = self.table.dumps().decode("utf-8").splitlines()
        self.assertEqual(lines[1], "\t".join(self.row))

    def test_load_without_probabilities(self):
        header = self.header[:7] + self.header[13:]
        row = self.row[:7] + self.row[13:]
        rows = ["\t".join(x) for x in (header, row)]
        data = "\n".join(rows).encode("utf-8")
        table = ClusterTable.loads(data)

    def test_dump_without_probabilities(self):
        row = self.row[:7] + self.row[13:]
        header = self.header[:7] + self.header[13:]
        rows = ["\t".join(x) for x in (header, row)]
        data = "\n".join(rows).encode("utf-8")
        table = ClusterTable.loads(data)
        lines = table.dumps().decode().splitlines()
        self.assertEqual(lines[1], "\t".join(self.row[:7] + self.row[13:]))
