"""Test `gecco.model.GeneTable` objects.
"""

import itertools
import io
import math
import os
import unittest
import warnings
from unittest import mock

import Bio.SeqIO
import polars
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gecco.model import Gene, GeneTable, Strand, Protein, Domain, _UnknownSeq


class TestGeneTable(unittest.TestCase):

    def test_from_genes(self):
        genes = [
            Gene(
                source=SeqRecord(id="seq1", seq=_UnknownSeq()),
                start=1,
                end=100,
                strand=Strand.Coding,
                protein=Protein(
                    id="seq1_1",
                    seq=None,
                    domains=[
                        Domain(
                            name="PF00106",
                            start=1,
                            end=20,
                            hmm="Pfam",
                            i_evalue=0.0,
                            pvalue=0.0,
                            probability=0.8
                        ),
                    ]
                )
            ),
            Gene(
                source=SeqRecord(id="seq2", seq=_UnknownSeq()),
                start=3,
                end=300,
                strand=Strand.Coding,
                protein=Protein(
                    id="seq2_1",
                    seq=None,
                    domains=[]
                ),
                _probability=0.6
            ),
            Gene(
                source=SeqRecord(id="seq2", seq=_UnknownSeq()),
                start=4,
                end=49,
                strand=Strand.Reverse,
                protein=Protein(
                    id="seq2_2",
                    seq=None,
                    domains=[],
                )
            ),
        ]

        table = GeneTable.from_genes(genes)
        self.assertEqual(len(table), 3)

        self.assertEqual(table.sequence_id[0], "seq1")
        self.assertEqual(table.protein_id[0], "seq1_1")
        self.assertEqual(table.start[0], 1)
        self.assertEqual(table.end[0], 100)
        self.assertEqual(table.strand[0], "+")
        self.assertEqual(table.average_p[0], 0.8)
        self.assertEqual(table.max_p[0], 0.8)

        self.assertEqual(table.sequence_id[1], "seq2")
        self.assertEqual(table.protein_id[1], "seq2_1")
        self.assertEqual(table.start[1], 3)
        self.assertEqual(table.end[1], 300)
        self.assertEqual(table.strand[1], "+")
        self.assertEqual(table.average_p[1], 0.6)
        self.assertEqual(table.max_p[1], 0.6)

        self.assertEqual(table.sequence_id[2], "seq2")
        self.assertEqual(table.protein_id[2], "seq2_2")
        self.assertEqual(table.start[2], 4)
        self.assertEqual(table.end[2], 49)
        self.assertEqual(table.strand[2], "-")
        self.assertTrue(math.isnan(table.average_p[2]))
        self.assertTrue(math.isnan(table.max_p[2]))

    def test_to_genes(self):
        table = GeneTable(polars.DataFrame(dict(
            sequence_id=["seq1", "seq2", "seq2"],
            protein_id=["seq1_1", "seq2_1", "seq2_2"],
            start=[100, 200, 300],
            end=[160, 260, 360],
            strand=["+", "-", "+"],
            average_p=[0.6, 0.2, None],
            max_p=[0.8, 0.2, None],
        )))

        genes = list(table.to_genes())

        self.assertEqual(genes[0].source.id, "seq1")
        self.assertEqual(genes[0].protein.id, "seq1_1")
        self.assertEqual(genes[0].start, 100)
        self.assertEqual(genes[0].end, 160)
        self.assertEqual(genes[0].strand, Strand.Coding)
        self.assertEqual(genes[0].average_probability, 0.6)

        self.assertEqual(genes[1].source.id, "seq2")
        self.assertEqual(genes[1].protein.id, "seq2_1")
        self.assertEqual(genes[1].start, 200)
        self.assertEqual(genes[1].end, 260)
        self.assertEqual(genes[1].strand, Strand.Reverse)
        self.assertEqual(genes[1].average_probability, 0.2)

        self.assertEqual(genes[2].source.id, "seq2")
        self.assertEqual(genes[2].protein.id, "seq2_2")
        self.assertEqual(genes[2].start, 300)
        self.assertEqual(genes[2].end, 360)
        self.assertEqual(genes[2].strand, Strand.Coding)
        self.assertIs(genes[2].average_probability, None)

    def test_dump(self):
        table = GeneTable(polars.DataFrame(dict(
            sequence_id=["seq_1", "seq_2", "seq_2"],
            protein_id=["seq1_1", "seq2_1", "seq2_2"],
            start=[100, 200, 300],
            end=[160, 260, 360],
            strand=["+", "-", "+"],
            average_p=[0.6, 0.2, math.nan],
            max_p=[0.8, 0.2, math.nan],
        )))

        lines = table.dumps().decode("utf-8").splitlines()

        self.assertEqual(
            lines[0],
            "\t".join(["sequence_id", "protein_id", "start", "end", "strand", "average_p", "max_p"])
        )
        self.assertEqual(
            lines[1],
            "\t".join(["seq_1", "seq1_1", "100", "160", "+", "0.6", "0.8"])
        )
        self.assertEqual(
            lines[2],
            "\t".join(["seq_2", "seq2_1", "200", "260", "-", "0.2", "0.2"])
        )
        self.assertEqual(
            lines[3],
            "\t".join(["seq_2", "seq2_2", "300", "360", "+", "", ""])
        )

    def test_dump_no_probability(self):
        table = GeneTable(polars.DataFrame(dict(
            sequence_id=["seq_1", "seq_2"],
            protein_id=["seq1_1", "seq2_1"],
            start=[100, 200],
            end=[160, 260],
            strand=["+", "-"],
            #average_p=[math.nan, math.nan],
            #max_p=[math.nan, math.nan],
        )))

        lines = table.dumps().decode("utf-8").splitlines()

        self.assertEqual(
            lines[0],
            "\t".join(["sequence_id", "protein_id", "start", "end", "strand"])
        )
        self.assertEqual(
            lines[1],
            "\t".join(["seq_1", "seq1_1", "100", "160", "+"])
        )
        self.assertEqual(
            lines[2],
            "\t".join(["seq_2", "seq2_1", "200", "260", "-"])
        )

    def test_load(self):
        lines = "\n".join([
            "\t".join(["sequence_id", "protein_id", "start", "end", "strand", "average_p", "max_p"]),
            "\t".join(["seq1", "seq1_1", "100", "160", "+", "0.6", "0.8"]),
            "\t".join(["seq2", "seq2_1", "200", "260", "-", "", ""]),
        ]).encode("utf-8")

        table = GeneTable.loads(lines)

        self.assertEqual(table.sequence_id[0], "seq1")
        self.assertEqual(table.protein_id[0], "seq1_1")
        self.assertEqual(table.start[0], 100)
        self.assertEqual(table.end[0], 160)
        self.assertEqual(table.strand[0], "+")
        self.assertEqual(table.average_p[0], 0.6)
        self.assertEqual(table.max_p[0], 0.8)

        self.assertEqual(table.sequence_id[1], "seq2")
        self.assertEqual(table.protein_id[1], "seq2_1")
        self.assertEqual(table.start[1], 200)
        self.assertEqual(table.end[1], 260)
        self.assertEqual(table.strand[1], "-")
        self.assertTrue(math.isnan(table.average_p[1]))
        self.assertTrue(math.isnan(table.max_p[1]))
