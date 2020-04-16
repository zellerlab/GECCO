"""Test `gecco.preprocessing` members.
"""

import unittest

import pandas
import gecco.crf.preprocessing


class TestSingleFeatureExtraction(unittest.TestCase):

    def test_single_column(self):
        extract = gecco.crf.preprocessing.extract_single_features
        df = pandas.DataFrame(
            columns=["protein_id", "domain", "p_domain", "BGC"],
            data=[
                ["prot1", "A", 0.5, "1"],
                ["prot1", "B", 1.0, "1"],
                ["prot2", "C", 0.8, "0"]
            ]
        )

        X, Y = extract(df, ["domain"], ["p_domain"], None)
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}])
        self.assertIs(Y, None)

        X, Y = extract(df, ["domain"], ["p_domain"], "BGC")
        self.assertEqual(X, [{"A": 0.5}, {"B": 1.0}, {"C": 0.8}])
        self.assertEqual(Y, ["1", "1", "0"])
