"""Test `gecco.crf.preprocessing` members.
"""

import unittest
import warnings

import pandas
import gecco.crf.preprocessing


class TestSingleFeatureExtraction(unittest.TestCase):

    # fmt: off
    def setUp(self):
        warnings.simplefilter("ignore")
        self.extract = gecco.crf.preprocessing.extract_single_features
        self.data = pandas.DataFrame(
            columns=["protein_id", "domain", "p_domain", "species", "p_species", "BGC"],
            data=[
                ["prot1", "A", 0.5, "NCBITaxon:1423", 0.95, "1"],
                ["prot1", "B", 1.0, "NCBITaxon:1423", 0.95, "1"],
                ["prot2", "C", 0.8, "NCBITaxon:1354003", 0.65, "0"],
            ]
        )

    def tearDown(self):
        warnings.simplefilter(warnings.defaultaction)

    def test_single_column_x_only(self):
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], None)
        self.assertIs(Y, None)
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}])

    def test_single_column(self):
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], "BGC")
        self.assertEqual(Y, ["1", "1", "0"])
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}])

    def test_multiple_columns_x_only(self):
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], None)
        self.assertIs(Y, None)
        self.assertEqual(X, [
            {"A": 0.5, "NCBITaxon:1423": 0.95},
            {"B": 1.0, "NCBITaxon:1423": 0.95},
            {"C": 0.8, "NCBITaxon:1354003": 0.65}
        ])

    def test_multiple_columns(self):
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], "BGC")
        self.assertEqual(Y, ["1", "1", "0"])
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}])


class TestGroupFeatureExtraction(unittest.TestCase):

    # fmt: off
    def setUp(self):
        warnings.simplefilter("ignore")
        self.extract = gecco.crf.preprocessing.extract_group_features
        self.data = pandas.DataFrame(
            columns=["seq_id", "protein_id", "domain", "p_domain", "species", "p_species", "BGC"],
            data=[
                ["seq1", "prot1", "A", 0.5, "NCBITaxon:1423", 0.95, "1"],
                ["seq1", "prot1", "B", 1.0, "NCBITaxon:1423", 0.95, "1"],
                ["seq1", "prot2", "C", 0.8, "NCBITaxon:1354003", 0.3, "0"],
                ["seq2", "prot3", "D", 0.4, "NCBITaxon:984897", 0.5, "0"],
            ]
        )

    def tearDown(self):
        warnings.simplefilter(warnings.defaultaction)

    def test_single_column_x_only(self):
        # group on protein id
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], "protein_id", None)
        self.assertIs(Y, None)
        self.assertEqual(X, [{"A": 0.5, "B": 1.0}, {"C": 0.8}, {"D": 0.4}])
        # group on sequence id
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], "seq_id", None)
        self.assertIs(Y, None)
        self.assertEqual(X, [{"A": 0.5, "B": 1.0, "C": 0.8}, {"D": 0.4}])

    def test_single_column(self):
        # group on protein id
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], "protein_id", "BGC")
        self.assertEqual(Y, ["1", "0", "0"])
        self.assertEqual(X, [{"A": 0.5, "B": 1.0}, {"C": 0.8}, {"D": 0.4}])
        # group on sequence id
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], "seq_id", "BGC")
        self.assertEqual(len(Y), 2)
        self.assertEqual(X, [{"A": 0.5, "B": 1.0, "C": 0.8}, {"D": 0.4}])

    def test_mixed_labels(self):
        # if mixed class labels withing a group, we should get a warning
        warnings.simplefilter("error")
        with self.assertRaises(UserWarning):
            self.extract(self.data, ["domain"], ["p_domain"], "seq_id", "BGC")
        with self.assertRaises(UserWarning):
            self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], "seq_id", "BGC")

    def test_multiple_columns_x_only(self):
        """Check that feature extraction from multiple column works.
        """
        # group on protein id
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], "protein_id", None)
        self.assertIs(Y, None)
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95},
            {"C": 0.8, "NCBITaxon:1354003": 0.3},
            {"D": 0.4, "NCBITaxon:984897": 0.5}
        ])
        # group on sequence id
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], "seq_id", None)
        self.assertIs(Y, None)
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95, "C": 0.8, "NCBITaxon:1354003": 0.3},
            {"D": 0.4, "NCBITaxon:984897": 0.5}
        ])

    def test_multiple_columns(self):
        # group on protein id
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], "protein_id", "BGC")
        self.assertEqual(Y, ["1", "0", "0"])
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95},
            {"C": 0.8, "NCBITaxon:1354003": 0.3},
            {"D": 0.4, "NCBITaxon:984897": 0.5}
        ])
        # group on sequence id
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], "seq_id", "BGC")
        self.assertEqual(len(Y), 2)
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95, "C": 0.8, "NCBITaxon:1354003": 0.3},
            {"D": 0.4, "NCBITaxon:984897": 0.5}
        ])


class TestOverlapFeatureExtraction(unittest.TestCase):

    # fmt: off
    def setUp(self):
        warnings.simplefilter("ignore")
        self.extract = gecco.crf.preprocessing.extract_overlapping_features
        self.data = pandas.DataFrame(
            columns=["protein_id", "domain", "p_domain", "species", "p_species", "BGC"],
            data=[
                ["prot1", "A", 0.5, "NCBITaxon:1423", 0.95, "1"],
                ["prot1", "B", 1.0, "NCBITaxon:1423", 0.95, "1"],
                ["prot2", "C", 0.8, "NCBITaxon:1354003", 0.3, "0"],
                ["prot3", "D", 0.4, "NCBITaxon:984897", 0.5, "0"],
            ]
        )

    def tearDown(self):
        warnings.simplefilter(warnings.defaultaction)

    def test_single_column_x_only(self):
        # overlap = 0
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], 0, None)
        self.assertIs(Y, None)
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}, {"D": 0.4}])
        # overlap = 1
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], 1, None)
        self.assertIs(Y, None)
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0},
            {"A": 0.5, "B": 1.0, "C": 0.8},
            {"B": 1.0, "C": 0.8, "D": 0.4},
            {"C": 0.8, "D": 0.4},
        ])
        # overlap = 2
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], 2, None)
        self.assertIs(Y, None)
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "C": 0.8},
            {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
            {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
            {"B": 1.0, "C": 0.8, "D": 0.4},
        ])

    def test_single_column(self):
        # overlap = 0
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], 0, "BGC")
        self.assertEqual(Y, ["1", "1", "0", "0"])
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}, {"D": 0.4}])
        # overlap = 1
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], 1, "BGC")
        self.assertEqual(Y, ["1", "1", "0", "0"])
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0},
            {"A": 0.5, "B": 1.0, "C": 0.8},
            {"B": 1.0, "C": 0.8, "D": 0.4},
            {"C": 0.8, "D": 0.4},
        ])
        # overlap = 2
        X, Y = self.extract(self.data, ["domain"], ["p_domain"], 2, "BGC")
        self.assertEqual(Y, ["1", "1", "0", "0"])
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "C": 0.8},
            {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
            {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
            {"B": 1.0, "C": 0.8, "D": 0.4},
        ])

    def test_multiple_columns_x_only(self):
            # overlap = 0
            X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 0, None)
            self.assertIs(Y, None)
            self.assertEqual(X, [
                {"A": 0.5, "NCBITaxon:1423": 0.95},
                {"B": 1.0, "NCBITaxon:1423": 0.95},
                {"C": 0.8, "NCBITaxon:1354003": 0.3},
                {"D": 0.4, "NCBITaxon:984897": 0.5}
            ])
            # overlap = 1
            X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 1, None)
            self.assertIs(Y, None)
            self.assertEqual(X, [
                {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95},
                {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
                {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
                {"C": 0.8, "D": 0.4, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
            ])
            # overlap = 2
            X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 2, None)
            self.assertIs(Y, None)
            self.assertEqual(X, [
                {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
                {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
                {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
                {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
            ])

    def test_multiple_columns(self):
        # overlap = 0
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 0, "BGC")
        self.assertEqual(Y, ["1", "1", "0", "0"])
        self.assertEqual(X, [
            {"A": 0.5, "NCBITaxon:1423": 0.95},
            {"B": 1.0, "NCBITaxon:1423": 0.95},
            {"C": 0.8, "NCBITaxon:1354003": 0.3},
            {"D": 0.4, "NCBITaxon:984897": 0.5}
        ])
        # overlap = 1
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 1, "BGC")
        self.assertEqual(Y, ["1", "1", "0", "0"])
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95},
            {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
            {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
            {"C": 0.8, "D": 0.4, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
        ])
        # overlap = 2
        X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 2, "BGC")
        self.assertEqual(Y, ["1", "1", "0", "0"])
        self.assertEqual(X, [
            {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
            {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
            {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
            {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
        ])
