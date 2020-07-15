import collections
import csv
import enum
import typing
from typing import Iterable, List, Optional, Sequence
from dataclasses import dataclass

import Bio.Alphabet
import numpy
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from . import __version__

if typing.TYPE_CHECKING:
    import io
    CsvWriter = type(csv.writer(io.StringIO()))


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
    i_evalue: float = 0.0


@dataclass
class Protein:
    id: str
    seq: Seq
    domains: Optional[List[Domain]] = None

    def to_record(self) -> SeqRecord:
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


@dataclass
class Cluster:
    id: str
    genes: List[Gene]

    type: str = "Other"
    type_probability: float = 0.0

    @property
    def seq_id(self):
        return self.genes[0].seq_id

    @property
    def start(self):
        return min(gene.start for gene in self.genes)

    @property
    def end(self):
        return max(gene.end for gene in self.genes)

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


def clusters_to_csv(clusters: Iterable["Cluster"], writer: "CsvWriter") -> int:
    # fmt: off
    writer.writerow([
        "sequence_id", "BGC_id", "start", "end", "average_p", "max_p",
        "BGC_type", "BGC_type_p", "proteins", "domains",
    ])

    written = 0
    for cluster in clusters:
        probs = numpy.array([ gene.probability for gene in cluster.genes ])
        writer.writerow([
            cluster.seq_id,
            cluster.id,
            cluster.start,
            cluster.end,
            probs.mean(),
            probs.max(),
            cluster.type,
            cluster.type_probability,
            ";".join([gene.id for gene in cluster.genes]),
            ";".join(sorted({
                domain.name
                for gene in cluster.genes
                for domain in gene.protein.domains
            })),
        ])
        written += 1
    return written
