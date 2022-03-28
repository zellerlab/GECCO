"""Test `gecco.model.GeneTable` objects.
"""

import itertools
import io
import os
import unittest
import warnings
from unittest import mock

import Bio.SeqIO
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

        rows = list(table)
        self.assertEqual(len(table), 3)
        self.assertEqual(len(rows), 3)

        row1 = rows[0]
        self.assertEqual(row1.sequence_id, "seq1")
        self.assertEqual(row1.protein_id, "seq1_1")
        self.assertEqual(row1.start, 1)
        self.assertEqual(row1.end, 100)
        self.assertEqual(row1.strand, "+")
        self.assertEqual(row1.average_p, 0.8)
        self.assertEqual(row1.max_p, 0.8)

        row2 = rows[1]
        self.assertEqual(row2.sequence_id, "seq2")
        self.assertEqual(row2.protein_id, "seq2_1")
        self.assertEqual(row2.start, 3)
        self.assertEqual(row2.end, 300)
        self.assertEqual(row2.strand, "+")
        self.assertEqual(row2.average_p, 0.6)
        self.assertEqual(row2.max_p, 0.6)

        row3 = rows[2]
        self.assertEqual(row3.sequence_id, "seq2")
        self.assertEqual(row3.protein_id, "seq2_2")
        self.assertEqual(row3.start, 4)
        self.assertEqual(row3.end, 49)
        self.assertEqual(row3.strand, "-")
        self.assertIs(row3.average_p, None)
        self.assertIs(row3.max_p, None)

    def test_to_genes(self):
        table = GeneTable(
            sequence_id=["seq1", "seq2", "seq2"],
            protein_id=["seq1_1", "seq2_1", "seq2_2"],
            start=[100, 200, 300],
            end=[160, 260, 360],
            strand=["+", "-", "+"],
            average_p=[0.6, 0.2, None],
            max_p=[0.8, 0.2, None],
        )

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
        table = GeneTable(
            sequence_id=["seq_1", "seq_2", "seq_2"],
            protein_id=["seq1_1", "seq2_1", "seq2_2"],
            start=[100, 200, 300],
            end=[160, 260, 360],
            strand=["+", "-", "+"],
            average_p=[0.6, 0.2, None],
            max_p=[0.8, 0.2, None],
        )

        buffer = io.StringIO()
        table.dump(buffer)
        lines = buffer.getvalue().splitlines()

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
        table = GeneTable(
            sequence_id=["seq_1", "seq_2"],
            protein_id=["seq1_1", "seq2_1"],
            start=[100, 200],
            end=[160, 260],
            strand=["+", "-"],
            average_p=[None, None],
            max_p=[None, None],
        )

        buffer = io.StringIO()
        table.dump(buffer)
        lines = buffer.getvalue().splitlines()

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
        ])

        table = GeneTable.load(io.StringIO(lines))

        self.assertEqual(table[0].sequence_id, "seq1")
        self.assertEqual(table[0].protein_id, "seq1_1")
        self.assertEqual(table[0].start, 100)
        self.assertEqual(table[0].end, 160)
        self.assertEqual(table[0].strand, "+")
        self.assertEqual(table[0].average_p, 0.6)
        self.assertEqual(table[0].max_p, 0.8)

        self.assertEqual(table[1].sequence_id, "seq2")
        self.assertEqual(table[1].protein_id, "seq2_1")
        self.assertEqual(table[1].start, 200)
        self.assertEqual(table[1].end, 260)
        self.assertEqual(table[1].strand, "-")
        self.assertEqual(table[1].average_p, None)
        self.assertEqual(table[1].max_p, None)
