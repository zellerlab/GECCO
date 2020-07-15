import collections
import csv
import enum
import typing
from typing import Iterable, List, Optional, Sequence
from dataclasses import dataclass

import Bio.Alphabet
import numpy
import pandas
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from . import __version__


class Strand(enum.IntEnum):
    Coding = 1
    Reverse = -1

    @property
    def sign(self):
        return "+" if self is Strand.Coding else "-"


@dataclass
class Domain:
    name: str
    start: int
    end: int
    hmm: str
    i_evalue: float = 0.0


@dataclass
class Protein:
    id: str
    seq: Seq
    domains: List[Domain]

    def __init__(self, id: str, seq: Seq, domains: Optional[List[Domain]] = None):
        self.id = id
        self.seq = seq
        self.domains = domains or list()

    def to_record(self) -> SeqRecord:
        # FIXME: add domains
        return SeqRecord(self.seq, id=self.id, name=self.id)


@dataclass
class Gene:
    seq_id: str
    start: int
    end: int
    strand: Strand
    protein: Protein
    probability: Optional[float] = None

    def __hash__(self):
        return hash(self.id)

    @property
    def id(self):
        return self.protein.id

    def to_feature_table(self) -> pandas.DataFrame:
        return pandas.DataFrame(
            data = [
                (
                    self.seq_id,
                    self.id,
                    self.start,
                    self.end,
                    self.strand.sign,
                    domain.name,
                    domain.hmm,
                    domain.i_evalue,
                    1 - domain.i_evalue,
                    domain.start,
                    domain.end,
                )
                for domain in self.protein.domains
            ],
            columns=[
                "sequence_id",
                "protein_id",
                "start",
                "end",
                "strand",
                "domain",
                "hmm",
                "i_Evalue",
                "rev_i_Evalue",
                "domain_start",
                "domain_end",
            ],
        )


@dataclass
class Cluster:
    id: str
    genes: List[Gene]

    type: str = "Other"
    type_probability: float = 0.0

    def __init__(
        self,
        id: str,
        genes: Optional[List[Gene]] = None,
        type: str = "Unknown",
        type_probability: float = 0.0
    ):
        self.id = id
        self.genes = genes or list()
        self.type = type
        self.type_probability = type_probability

    @property
    def seq_id(self):
        return self.genes[0].seq_id

    @property
    def start(self):
        return min(gene.start for gene in self.genes)

    @property
    def end(self):
        return max(gene.end for gene in self.genes)

    @property
    def average_probability(self):
        p = [g.probability for g in self.genes if g.probability is not None]
        return sum(p) / len(p)

    @property
    def maximum_probability(self):
        p = [g.probability for g in self.genes if g.probability is not None]
        return max(p)

    # ---

    def domain_composition(self, all_possible: Optional[Sequence[str]] = None) -> numpy.ndarray:
        domains = [d for gene in self.genes for d in gene.protein.domains]
        names = numpy.array([domain.name for domain in domains])
        weights = numpy.array([1 - domain.i_evalue for domain in domains])
        counts = collections.Counter(names)
        if all_possible is None:
            all_possible = numpy.unique(names)
        composition = numpy.zeros(len(all_possible))
        for i, dom in enumerate(all_possible):
            n = counts[dom]
            w = weights[names == dom].mean() if n > 0 else 0
            composition[i] = n * w
        return composition / (composition.sum() or 1)

    # ---

    def to_record(self, sequence: SeqRecord) -> SeqRecord:
        if sequence.id != self.genes[0].seq_id:
            raise ValueError(
                f"invalid sequence {sequence.id!r} "
                f"(expected {self.genes[0].seq_id!r})"
            )

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
            cds = SeqFeature(location=loc, type="cds")
            cds.qualifiers["label"] = gene.id
            cds.qualifiers["inference"] = ["EXISTENCE:ab initio prediction:PRODIGAL:2.6"]
            cds.qualifiers["translation"] = str(gene.protein.seq)
            bgc.features.append(cds)

            # write domains as /misc_feature annotations
            for domain in gene.protein.domains:
                dom_start = start + domain.start * 3
                dom_end = start + domain.end * 3
                loc = FeatureLocation(dom_start, dom_end, strand=int(gene.strand))
                misc = SeqFeature(location=loc, type="misc_feature")
                misc.qualifiers["label"] = domain.name
                bgc.features.append(misc)

        return bgc

    def to_cluster_table(self) -> pandas.DataFrame:
        return pandas.DataFrame(
            data = [
                (
                    self.seq_id,
                    self.id,
                    self.start,
                    self.end,
                    self.average_probability,
                    self.maximum_probability,
                    self.type,
                    self.type_probability,
                    ";".join([gene.id for gene in self.genes]),
                    ";".join(sorted({
                        domain.name
                        for gene in self.genes
                        for domain in gene.protein.domains
                    })),
                )
            ],
            columns=[
                "sequence_id", "BGC_id", "start", "end", "average_p", "max_p",
                "BGC_type", "BGC_type_p", "proteins", "domains"
            ],
        )
