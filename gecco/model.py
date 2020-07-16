"""Data layer classes storing information needed for BGC detection.
"""

import csv
import enum
import re
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


class Hmm(typing.NamedTuple):

    id: str
    version: str
    url: str
    path: str
    relabel_with: Optional[str] = None

    def relabel(self, domain: str) -> str:
        if self.relabel_with is None:
            return domain
        before, after = re.match("^s/(.*)/(.*)/$", self.relabel_with).groups()  # type: ignore
        regex = re.compile(before)
        return regex.sub(after, domain)


class Strand(enum.IntEnum):
    """A flag to declare on which DNA strand a gene is located.
    """

    Coding = 1
    Reverse = -1

    @property
    def sign(self) -> str:
        """`str`: The strand as a single sign (``+`` or ``-``).
        """
        return "+" if self is Strand.Coding else "-"


class Domain(typing.NamedTuple):
    """A conserved region within a protein.

    Attributes:
        name (`str`): The accession of the protein domain in the source HMM.
        start (`int`): The start coordinate of the domain within the protein
            sequence (first amino-acid at 1).
        end (`int`): The end coordinate of the domain within the protein
            sequence (inclusive).
        hmm (`str`): The name of the HMM library this domain belongs to
            (e.g. ``Pfam``, ``Panther``).
        i_evalue (`float`): The independent e-value reported by ``hmmsearch``
            that measures how reliable the domain annotation is.

    """

    name: str
    start: int
    end: int
    hmm: str
    i_evalue: float = 0.0


@dataclass
class Protein:
    """A sequence of amino-acids translated from a gene.

    Attributes:
        id (`str`): The identifier of the protein.
        seq (`~Bio.Seq.Seq`): The sequence of amino-acids of this protein.
        domains (`list` of `~gecco.model.Domain`): A list of domains found
            in the protein sequence.

    """

    id: str
    seq: Seq
    domains: List[Domain]

    def __init__(self, id: str, seq: Seq, domains: Optional[List[Domain]] = None):  # noqa: D107
        self.id = id
        self.seq = seq
        self.domains = domains or list()

    def to_record(self) -> SeqRecord:
        """Convert the protein to a single record.
        """
        # FIXME: add domains
        return SeqRecord(self.seq, id=self.id, name=self.id)


@dataclass
class Gene:
    """A nucleotide sequence coding a protein.

    Attributes:
        source (`~Bio.SeqRecord.SeqRecord`): The DNA sequence this gene was
            found in, as a Biopython record.
        start (`int`): The index of the leftmost nucleotide of the gene within
            the source sequence, independent of the strandedness.
        end (`int`): The index of the rightmost nucleotide of the gene within
            the source sequence.
        protein (`~gecco.model.Protein`): The protein translated from this
            gene.
        probability (`float`, optional): The probability with which this
            protein is part of a biosynthetic gene cluster.

    """

    source: SeqRecord
    start: int
    end: int
    strand: Strand
    protein: Protein
    probability: Optional[float] = None

    @property
    def id(self) -> str:
        """`str`: The identifier of the gene (same as the protein identifier).
        """
        return self.protein.id

    def to_feature_table(self) -> pandas.DataFrame:
        """Convert this gene to a feature table listing domain annotations.

        The returned objects can be concatenated together to obtain a feature
        table that can be passed to a `~gecco.crf.ClusterCRF` instance.

        Returns:
            `pandas.DataFrame`: A dataframe listing all domain annotation
            in no particular order.

        """
        return pandas.DataFrame(
            data = [
                (
                    self.source.id,
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
    """A sequence of contiguous genes with biosynthetic activity.

    Attributes:
        id (`str`): The identifier of the gene cluster.
        genes (`list` of `~gecco.model.Gene`): A list of the genes belonging
            to this gene cluster.
        type (`str`): The putative type of product synthesized by this gene
            cluster, according to similarity with existing clusters.
        type_probability (`float`): The probability with which the BGC type
            was identified.

    """

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
    ):  # noqa: D107
        self.id = id
        self.genes = genes or list()
        self.type = type
        self.type_probability = type_probability

    @property
    def source(self) -> SeqRecord:  # type: ignore
        """`~Bio.SeqRecord.SeqRecord`: The sequence this cluster was found in.
        """
        return self.genes[0].source

    @property
    def start(self) -> int:
        """`int`: The start of this cluster in the source sequence.
        """
        return min(gene.start for gene in self.genes)

    @property
    def end(self) -> int:
        """`int`: The end of this cluster in the source sequence.
        """
        return max(gene.end for gene in self.genes)

    @property
    def average_probability(self) -> float:
        """`float`: The average of proteins probability of being biosynthetic.
        """
        p = [g.probability for g in self.genes if g.probability is not None]
        return sum(p) / len(p)

    @property
    def maximum_probability(self) -> float:
        """`float`: The highest of proteins probability of being biosynthetic.
        """
        p = [g.probability for g in self.genes if g.probability is not None]
        return max(p)

    # ---

    def domain_composition(self, all_possible: Optional[Sequence[str]] = None) -> numpy.ndarray:
        """Compute weighted domain composition with respect to ``all_possible``.

        Arguments:
            all_possible (sequence of `str`, optional): A sequence containing
                all domain names to consider when computing domain composition
                for the BGC. If `None` given, then only domains within the
                cluster are taken into account.

        Returns:
            `~numpy.ndarray`: A numerical array containing the relative domain
            composition of the BGC.

        """
        domains = [d for gene in self.genes for d in gene.protein.domains]
        names = numpy.array([domain.name for domain in domains])
        weights = numpy.array([1 - domain.i_evalue for domain in domains])
        if all_possible is None:
            all_possible = numpy.unique(names)
        composition = numpy.zeros(len(all_possible))
        for i, dom in enumerate(all_possible):
            composition[i] = numpy.sum(weights[names == dom])
        return composition / (composition.sum() or 1) # type: ignore

    # ---

    def to_record(self) -> SeqRecord:
        """Convert the cluster to a single record.

        Annotations of the source sequence are kept intact if they don't
        overlap with the cluster boundaries. Component genes are added on the
        record as *CDS* features. Annotated protein domains are added as
        *misc_feature*.

        """
        bgc = self.source[self.start:self.end]
        bgc.id = bgc.name = self.id
        bgc.seq.alphabet = Bio.Alphabet.generic_dna

        # copy sequence annotations
        bgc.annotations["topology"] = "linear"
        bgc.annotations["organism"] = self.source.annotations.get("organism")
        bgc.annotations["source"] = self.source.annotations.get("source")
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
        """Convert this cluster to a cluster table with a single row.

        The obtained `pandas.DataFrame` can be concatenated together with
        other tables, which allows to easily convert a sequence of
        `~gecco.model.Cluster` into a single `~pandas.DataFrame` quite easily.

        """
        return pandas.DataFrame(
            data = [
                (
                    self.source.id,
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

    def to_feature_table(self) -> pandas.DataFrame:
        """Convert this cluster to a feature table listing domain annotations.
        """
        return pandas.concat(map(Gene.to_feature_table, self.genes))
