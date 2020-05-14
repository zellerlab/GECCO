"""Test `gecco.crf` members.
"""

import itertools
import os
import unittest
import warnings
from unittest import mock

import Bio.SeqIO
import pandas
from gecco.crf import ClusterCRF
from gecco.crf.select import fisher_significance


class TestClusterCRF(unittest.TestCase):

    def test_fisher_selection_bug(self):
        """Check that marginals are emitted for rows without any selected feature.
        """
        # we have a dummy sequence with two proteins: protA and protB,
        # which have been annotated with domains A and B, resp. A C and D.
        crf = ClusterCRF("domain", "p_domain")#, feature_type="group")
        data = [pandas.DataFrame(
            columns=["sequence_id", "protein_id", "domain", "p_domain", "BGC"],
            data=[
                ("seqA", "protA", "A", 0.9, 1),
                ("seqA", "protA", "B", 0.7, 1),
                ("seqA", "protB", "A", 0.9, 0),
                ("seqA", "protB", "C", 0.6, 0),
                ("seqA", "protB", "D", 0.9, 0),
            ]
        )]

        # we patch `fisher_significance` to make sure ClusterCRF would favour
        # B, C, and D, since A is the least significant of the domains, but
        # because of the small number of samples a Fisher Exact Test would
        # always accept the null-hypothesis.
        with mock.patch("gecco.crf.fisher_significance") as _patch:
            _patch.return_value = {"A":1, "B":0, "C": 0, "D": 0}
            crf.fit(data, select=0.75)
            self.assertEqual(crf.significant_features["domain"], {"B", "C", "D"})

        # now, trying to predict with C and D should give us a non-BGC
        pred_data = [pandas.DataFrame(
            columns=["sequence_id", "protein_id", "domain", "p_domain"],
            data=[
                ("seqA", "protB", "C", 0.8),
                ("seqA", "protB", "D", 0.8),
            ]
        )]
        marginals = crf.predict_marginals(pred_data)
        self.assertEqual(len(marginals), 2)
        self.assertLessEqual(marginals.p_pred.mean(), 0.5)

        # and trying to predict with B should give us a BGC
        pred_data = [pandas.DataFrame(
            columns=["sequence_id", "protein_id", "domain", "p_domain"],
            data=[("seqA", "protB", "B", 0.8)]
        )]
        marginals = crf.predict_marginals(pred_data)
        self.assertEqual(len(marginals), 1)
        self.assertGreaterEqual(marginals.p_pred.mean(), 0.5)

        # but trying to predict with only A should give a prediction, although
        # the CRF has no clue what to do in this situation since it didn't
        # take into account samples annotated with A, so it will default to
        # a random prediction
        pred_data = [pandas.DataFrame(
            columns=["sequence_id", "protein_id", "domain", "p_domain"],
            data=[("seqA", "protB", "A", 0.7)]
        )]
        marginals = crf.predict_marginals(pred_data)
        self.assertEqual(len(marginals), 1)
        self.assertEqual(marginals.p_pred.mean(), 0.5)
