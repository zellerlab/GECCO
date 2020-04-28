"""Test `gecco.bgc` members.
"""

import io
import itertools
import os
import unittest
import warnings
from unittest import mock

import Bio.SeqIO
import pandas
import numpy
import gecco.bgc
from gecco.bgc import BGC, Protein


class TestProtein(unittest.TestCase):

    def test_coordinates(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6)
        self.assertEqual(p1.start, 100)
        self.assertEqual(p1.end, 200)

        p2 = Protein("seq1", 500, 200, "prot2", 0.8)
        self.assertEqual(p2.start, 200)
        self.assertEqual(p2.end, 500)

    def test_invalid_weight_array(self):
        with self.assertRaises(ValueError):
            Protein("seq1", 100, 200, "prot1", domains=["A", "B", "C"], weights=[0.5])

    def test_is_potential_cluster(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6)
        self.assertTrue(p1.is_potential_cluster())
        self.assertFalse(p1.is_potential_cluster(0.8))

        p2 = Protein("seq1", 100, 200, "prot1", 0.2)
        self.assertFalse(p2.is_potential_cluster())
        self.assertTrue(p2.is_potential_cluster(0.1))


class TestBGC(unittest.TestCase):

    def test_seq_id(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6)
        p2 = Protein("seq1", 500, 300, "prot2", 0.8)
        bgc = BGC([p1, p2])
        self.assertEqual(bgc.seq_id, "seq1")

    def test_coordinates(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6)
        p2 = Protein("seq1", 500, 300, "prot2", 0.8)

        bgc = BGC([p1, p2])
        self.assertEqual(bgc.start, 100)
        self.assertEqual(bgc.end, 500)
        bgc = BGC([p2, p1])
        self.assertEqual(bgc.start, 100)
        self.assertEqual(bgc.end, 500)

    def test_empty_bgc(self):
        self.assertRaises(ValueError, BGC, proteins=[])

    def test_is_valid_gecco(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6)
        p2 = Protein("seq1", 500, 200, "prot2", 0.8)
        bgc = BGC([p1, p2])
        self.assertTrue(bgc.is_valid(criterion="gecco"))

    def test_is_valid_antismash(self):
        bio_pfams = list(gecco.bgc.BIO_PFAMS)

        p1 = Protein("seq1", 10, 100, "prot1", 0.8, domains=bio_pfams[:2])
        p2 = Protein("seq1", 200, 1000, "prot2", 0.9,  domains=bio_pfams[:4])
        p3 = Protein("seq1", 800, 340, "prot3", 0.84,  domains=bio_pfams[10:15])
        p4 = Protein("seq1", 2000, 2200, "prot4", 0.8, domains=bio_pfams[20:21])
        p5 = Protein("seq1", 2300, 2360, "prot5", 0.82, domains=bio_pfams[24:28])
        p6 = Protein("seq1", 3000, 3500, "prot6", 0.81, domains=bio_pfams[21:23])
        bgc = BGC([p1, p2, p3, p4, p5, p6])

        # actually valid
        self.assertTrue(bgc.is_valid(criterion="antismash"))
        # not valid after we raised the threshold
        self.assertFalse(bgc._antismash_check(p_thresh=0.9))
        # not valid after we raise the bio Pfam count
        self.assertFalse(bgc._antismash_check(n_biopfams=20))
        # not valid after we raise the CDS count
        self.assertFalse(bgc._antismash_check(n_cds=bgc.prot_ids.size+1))

    def test_is_valid_unknown_criterion(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6, ["A", "B"])
        p2 = Protein("seq1", 500, 200, "prot2", 0.8, ["A", "C"])
        bgc = BGC([p1, p2], name="bgc1")
        self.assertRaises(ValueError, bgc.is_valid, criterion="unknown")

    def test_write_to_file(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6, ["A", "B"])
        p2 = Protein("seq1", 500, 200, "prot2", 0.8, ["A", "C"])
        bgc = BGC([p1, p2], name="bgc1")

        with io.StringIO() as buffer:
            bgc.write_to_file(buffer)
            self.assertEqual(
                buffer.getvalue(),
                "seq1\tbgc1\t100\t500\t0.7\t0.8\t\t\r\n"
            )

        with io.StringIO() as buffer:
            bgc.write_to_file(buffer, long=True)
            self.assertEqual(
                buffer.getvalue(),
                "seq1\tbgc1\t100\t500\t0.7\t0.8\t\t\tprot1;prot2\tA;B;A;C\r\n"
            )

    def test_domain_counts(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6, ["A", "B"])
        p2 = Protein("seq1", 500, 200, "prot2", 0.8, ["A", "C"])
        bgc = BGC([p1, p2], name="bgc1")

        self.assertEqual(list(bgc.domain_counts()), [2, 1, 1])
        self.assertEqual(list(bgc.domain_counts(["B", "C"])), [1, 1])
        self.assertEqual(list(bgc.domain_counts(["B", "A"])), [1, 2])

    def test_domain_composition(self):
        p1 = Protein("seq1", 100, 200, "prot1", 0.6, ["A", "B"])
        p2 = Protein("seq1", 500, 200, "prot2", 0.8, ["A", "C"])
        bgc = BGC([p1, p2], name="bgc1")

        self.assertEqual(list(bgc.domain_composition()), [0.5, 0.25, 0.25])
        self.assertEqual(list(bgc.domain_composition(["B", "C"])), [0.5, 0.5])
        self.assertEqual(list(bgc.domain_composition(["B", "C", "A"])), [0.25, 0.25, 0.5])
