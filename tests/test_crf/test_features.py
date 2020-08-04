"""Test `gecco.crf.preprocessing` members.
"""

import unittest
import warnings

import pandas
import gecco.crf.features
from gecco.model import Domain, Gene, Protein, Strand


class _TestFeatureExtraction(object):

    def setUp(self):
        warnings.simplefilter("ignore")
        self.genes = [
            Gene(
                source=None,
                start=0,
                end=1,
                strand=Strand.Coding,
                protein=Protein(
                    id="prot1",
                    seq=None,
                    domains=[
                        Domain("A", 0, 1, "test", 0.5, probability=1, qualifiers={}),
                        Domain("B", 0, 1, "test", 0.0, probability=1, qualifiers={})
                    ]
                )
            ),
            Gene(
                source=None,
                start=1,
                end=2,
                strand=Strand.Coding,
                protein=Protein(
                    id="prot2",
                    seq=None,
                    domains=[
                        Domain("C", 0, 1, "test", 0.2, probability=0, qualifiers={})
                    ]
                )
            )
        ]

    def tearDown(self):
        warnings.simplefilter(warnings.defaultaction)


class TestSingleFeatureExtraction(_TestFeatureExtraction, unittest.TestCase):

    extract = staticmethod(gecco.crf.features.extract_features_single)

    def test_extraction(self):
        X = self.extract(self.genes)
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}])


class TestGroupFeatureExtraction(_TestFeatureExtraction, unittest.TestCase):

    extract = staticmethod(gecco.crf.features.extract_features_group)

    def test_extraction(self):
        X = self.extract(self.genes)
        self.assertEqual(X, [{"A": 0.5, "B": 1.0}, {"C": 0.8}])


# class TestOverlapFeatureExtraction(unittest.TestCase):
#
#     # fmt: off
#     def setUp(self):
#         warnings.simplefilter("ignore")
#         self.extract = gecco.crf.preprocessing.extract_overlapping_features
#         self.data = pandas.DataFrame(
#             columns=["protein_id", "domain", "p_domain", "species", "p_species", "BGC"],
#             data=[
#                 ["prot1", "A", 0.5, "NCBITaxon:1423", 0.95, "1"],
#                 ["prot1", "B", 1.0, "NCBITaxon:1423", 0.95, "1"],
#                 ["prot2", "C", 0.8, "NCBITaxon:1354003", 0.3, "0"],
#                 ["prot3", "D", 0.4, "NCBITaxon:984897", 0.5, "0"],
#             ]
#         )
#
#     def tearDown(self):
#         warnings.simplefilter(warnings.defaultaction)
#
#     def test_single_column_x_only(self):
#         # overlap = 0
#         X, Y = self.extract(self.data, ["domain"], ["p_domain"], 0, None)
#         self.assertIs(Y, None)
#         self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}, {"D": 0.4}])
#         # overlap = 1
#         X, Y = self.extract(self.data, ["domain"], ["p_domain"], 1, None)
#         self.assertIs(Y, None)
#         self.assertEqual(X, [
#             {"A": 0.5, "B": 1.0},
#             {"A": 0.5, "B": 1.0, "C": 0.8},
#             {"B": 1.0, "C": 0.8, "D": 0.4},
#             {"C": 0.8, "D": 0.4},
#         ])
#         # overlap = 2
#         X, Y = self.extract(self.data, ["domain"], ["p_domain"], 2, None)
#         self.assertIs(Y, None)
#         self.assertEqual(X, [
#             {"A": 0.5, "B": 1.0, "C": 0.8},
#             {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
#             {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
#             {"B": 1.0, "C": 0.8, "D": 0.4},
#         ])
#
#     def test_single_column(self):
#         # overlap = 0
#         X, Y = self.extract(self.data, ["domain"], ["p_domain"], 0, "BGC")
#         self.assertEqual(Y, ["1", "1", "0", "0"])
#         self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}, {"D": 0.4}])
#         # overlap = 1
#         X, Y = self.extract(self.data, ["domain"], ["p_domain"], 1, "BGC")
#         self.assertEqual(Y, ["1", "1", "0", "0"])
#         self.assertEqual(X, [
#             {"A": 0.5, "B": 1.0},
#             {"A": 0.5, "B": 1.0, "C": 0.8},
#             {"B": 1.0, "C": 0.8, "D": 0.4},
#             {"C": 0.8, "D": 0.4},
#         ])
#         # overlap = 2
#         X, Y = self.extract(self.data, ["domain"], ["p_domain"], 2, "BGC")
#         self.assertEqual(Y, ["1", "1", "0", "0"])
#         self.assertEqual(X, [
#             {"A": 0.5, "B": 1.0, "C": 0.8},
#             {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
#             {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4},
#             {"B": 1.0, "C": 0.8, "D": 0.4},
#         ])
#
#     def test_multiple_columns_x_only(self):
#             # overlap = 0
#             X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 0, None)
#             self.assertIs(Y, None)
#             self.assertEqual(X, [
#                 {"A": 0.5, "NCBITaxon:1423": 0.95},
#                 {"B": 1.0, "NCBITaxon:1423": 0.95},
#                 {"C": 0.8, "NCBITaxon:1354003": 0.3},
#                 {"D": 0.4, "NCBITaxon:984897": 0.5}
#             ])
#             # overlap = 1
#             X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 1, None)
#             self.assertIs(Y, None)
#             self.assertEqual(X, [
#                 {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95},
#                 {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
#                 {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#                 {"C": 0.8, "D": 0.4, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#             ])
#             # overlap = 2
#             X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 2, None)
#             self.assertIs(Y, None)
#             self.assertEqual(X, [
#                 {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
#                 {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#                 {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#                 {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#             ])
#
#     def test_multiple_columns(self):
#         # overlap = 0
#         X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 0, "BGC")
#         self.assertEqual(Y, ["1", "1", "0", "0"])
#         self.assertEqual(X, [
#             {"A": 0.5, "NCBITaxon:1423": 0.95},
#             {"B": 1.0, "NCBITaxon:1423": 0.95},
#             {"C": 0.8, "NCBITaxon:1354003": 0.3},
#             {"D": 0.4, "NCBITaxon:984897": 0.5}
#         ])
#         # overlap = 1
#         X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 1, "BGC")
#         self.assertEqual(Y, ["1", "1", "0", "0"])
#         self.assertEqual(X, [
#             {"A": 0.5, "B": 1.0, "NCBITaxon:1423": 0.95},
#             {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
#             {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#             {"C": 0.8, "D": 0.4, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#         ])
#         # overlap = 2
#         X, Y = self.extract(self.data, ["domain", "species"], ["p_domain", "p_species"], 2, "BGC")
#         self.assertEqual(Y, ["1", "1", "0", "0"])
#         self.assertEqual(X, [
#             {"A": 0.5, "B": 1.0, "C": 0.8, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3},
#             {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#             {"A": 0.5, "B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#             {"B": 1.0, "C": 0.8, "D": 0.4, "NCBITaxon:1423": 0.95, "NCBITaxon:1354003": 0.3, "NCBITaxon:984897": 0.5},
#         ])
