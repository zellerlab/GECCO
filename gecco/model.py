import enum
import typing
from typing import List, Optional

import Bio.Alphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from . import __version__

# ----------------------------------------------------------------------------

class Strand(enum.IntEnum):
    Coding = 1
    Reverse = -1

    @property
    def sign(self):
        return "+" if self is Strand.Coding else "-"


# ----------------------------------------------------------------------------


class Protein(typing.NamedTuple):
    id: str
    seq: Seq

    def to_record(self) -> SeqRecord:
        return SeqRecord(self.seq, id=self.id, name=self.id)


class Gene(typing.NamedTuple):
    seq_id: str
    start: int
    end: int
    strand: Strand
    protein: Protein

    @property
    def id(self):
        return self.protein.id


# ----------------------------------------------------------------------------


class Domain(typing.NamedTuple):
    name: str
    start: int
    end: int
    i_evalue: float


class ClusterProtein(typing.NamedTuple, Protein):
    id: str
    seq: Seq
    domains: List[Domain]


class ClusterGene(typing.NamedTuple, Gene):
    seq_id: str
    start: int
    end: int
    strand: Strand
    protein: ClusterProtein

    @property
    def id(self):
        return self.protein.id


class Cluster(typing.NamedTuple):
    id: str
    genes: List[ClusterGene]

    type: str = "Other"
    type_prob: float = 0

    @property
    def seq_id(self):
        return self.genes[0].seq_id

    @property
    def start(self):
        return min(gene.start for gene in self.genes)

    @property
    def end(self):
        return max(gene.end for gene in self.genes)

    def is_valid(self, criterion: str = "gecco"):
        return True # TODO

    def to_record(self, sequence: SeqRecord):
        bgc = sequence[self.start:self.end]
        bgc.id = bgc.name = self.id
        bgc.seq.alphabet = Bio.Alphabet.generic_dna

        # copy sequence annotations
        bgc.annotations["topology"] = "linear"
        bgc.annotations["organism"] = sequence.annotations.get("organism")
        bgc.annotations["source"] = sequence.annotations.get("source")
        bgc.annotations["comment"] = ["Detected with GECCO v{}".format(__version__)]

        # add proteins as CDS features
        for gene in self.genes:
            # translate gene location
            start = gene.start - self.start - 1
            end = gene.end - self.start
            loc = FeatureLocation(start=start, end=end, strand=int(gene.strand))
            # write protein as a /cds GenBank record
            cds = SeqFeature(location=loc, type="cds", id=gene.id)
            cds.qualifiers["label"] = gene.id
            cds.qualifiers["inference"] = ["EXISTENCE:ab initio prediction:PRODIGAL:2.6"]
            cds.qualifiers["translation"] = str(gene.protein.seq)
            bgc.features.append(cds)

            # TODO: include domain annotations

        return bgc
