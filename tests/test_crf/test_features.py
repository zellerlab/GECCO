"""Test `gecco.crf.preprocessing` members.
"""

import unittest
import warnings

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
