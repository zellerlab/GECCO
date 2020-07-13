import enum
import typing
from typing import List, Optional

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ----------------------------------------------------------------------------

class Strand(enum.Enum):
    Coding = enum.auto()
    Reverse = enum.auto()

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

class ClusterProtein(Protein):
    pass

class ClusterGene(Gene):
    protein: ClusterProtein
