"""Data layer classes storing information needed for BGC detection.
"""

import csv
import enum
import itertools
import operator
import re
import typing
from collections.abc import Sized
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, TextIO, NamedTuple, Union, Iterator

import Bio.Alphabet
import numpy
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from . import __version__
from ._base import Dumpable


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


@dataclass(frozen=True)
class Domain:
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
        probability (`float`, optional): The probability that this domain
            is part of a BGC, or `None` if no prediction has been made yet.
        qualifiers (`dict`, optional): A dictionary of feature qualifiers that
            is added to the `~Bio.SeqFeature.SeqFeature` built from this
            `Domain`.

    """

    name: str
    start: int
    end: int
    hmm: str
    i_evalue: float
    probability: Optional[float] = None
    qualifiers: Mapping[str, Union[str, List[str]]] = field(default_factory=dict)

    def with_probability(self, probability: Optional[float]) -> "Domain":
        """Copy the current domain and assign it a BGC probability.
        """
        return Domain(
            self.name, self.start, self.end, self.hmm, self.i_evalue,
            probability, self.qualifiers
        )

    def to_seq_feature(self, protein_coordinates: bool = False) -> SeqFeature:
        """Convert the domain to a single feature.

        Arguments:
            protein_coordinates (`bool`): Set to `True` for the feature
                coordinates to be given in amino-acids, or to `False` in
                nucleotides.

        """
        stride = 1 if protein_coordinates else 3
        loc = FeatureLocation(0, (self.end - self.start)*stride)
        qualifiers = dict(self.qualifiers)
        qualifiers.setdefault("standard_name", self.name)
        return SeqFeature(location=loc, type="misc_feature", qualifiers=qualifiers)


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
    domains: List[Domain] = field(default_factory=list)

    def to_seq_record(self) -> SeqRecord:
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
        qualifiers (`dict`, optional): A dictionary of feature qualifiers that
            is added to the `~Bio.SeqFeature.SeqFeature` built from this
            `Gene`.

    """

    source: SeqRecord
    start: int
    end: int
    strand: Strand
    protein: Protein
    qualifiers: Mapping[str, Union[str, List[str]]] = field(default_factory=dict)

    @property
    def id(self) -> str:
        """`str`: The identifier of the gene (same as the protein identifier).
        """
        return self.protein.id

    @property
    def average_probability(self) -> Optional[float]:
        """`float`: The average of domain probabilities of being biosynthetic.
        """
        p = [d.probability for d in self.protein.domains if d.probability is not None]
        return sum(p) / len(p) if p else None

    @property
    def maximum_probability(self) -> Optional[float]:
        """`float`: The highest of domain probabilities of being biosynthetic.
        """
        p = [d.probability for d in self.protein.domains if d.probability is not None]
        return max(p) if p else None

    # ---

    def to_seq_feature(self) -> SeqFeature:
        """Convert the gene to a single feature.
        """
        # NB(@althonos): we use inclusive 1-based ranges in the data model
        # but Biopython expects 0-based ranges with exclusive ends
        end = self.end - self.start + 1
        loc = FeatureLocation(start=0, end=end, strand=int(self.strand))
        qualifiers = dict(self.qualifiers)
        qualifiers.setdefault("locus_tag", self.protein.id)
        qualifiers.setdefault("translation", str(self.protein.seq))
        return SeqFeature(location=loc, type="cds", qualifiers=qualifiers)


@dataclass
class Cluster:
    """A sequence of contiguous genes with biosynthetic activity.

    Attributes:
        id (`str`): The identifier of the gene cluster.
        genes (`list` of `~gecco.model.Gene`): A list of the genes belonging
            to this gene cluster.
        types (`list` of `str`): The list of the putative types of product
            synthesized by this gene cluster, according to similarity in
            domain composition with curated clusters.
        types_probabilities (`list` of `float`): The probability with which
            each BGC type was identified (same dimension as the ``types``
            attribute).

    """

    id: str
    genes: List[Gene]
    types: List[str]
    types_probabilities: List[float]

    def __init__(
        self,
        id: str,
        genes: Optional[List[Gene]] = None,
        types: Optional[List[str]] = None,
        types_probabilities: Optional[List[float]] = None,
    ):  # noqa: D107
        self.id = id
        self.genes = genes or list()
        self.types = types or list()
        self.types_probabilities = types_probabilities or list()

        if len(self.types) != len(self.types_probabilities):
            err = "type and type probability lists must have the same dimensions"
            raise ValueError(err)


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
    def average_probability(self) -> Optional[float]:
        """`float`: The average of proteins probability of being biosynthetic.
        """
        p = [g.average_probability for g in self.genes if g.average_probability is not None]
        return sum(p) / len(p) if p else None

    @property
    def maximum_probability(self) -> Optional[float]:
        """`float`: The highest of proteins probability of being biosynthetic.
        """
        p = [g.maximum_probability for g in self.genes if g.maximum_probability is not None]
        return max(p) if p else None

    # ---

    def domain_composition(
        self, all_possible: Optional[Sequence[str]] = None
    ) -> numpy.ndarray:
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
        unique_names = set(names)
        if all_possible is None:
            all_possible = numpy.unique(names)
        composition = numpy.zeros(len(all_possible))
        for i, dom in enumerate(all_possible):
            if dom in unique_names:
                composition[i] = numpy.sum(weights[names == dom])
        return composition / (composition.sum() or 1)  # type: ignore

    # ---

    def to_seq_record(self) -> SeqRecord:
        """Convert the cluster to a single record.

        Annotations of the source sequence are kept intact if they don't
        overlap with the cluster boundaries. Component genes are added on the
        record as *CDS* features. Annotated protein domains are added as
        *misc_feature*.

        """
        # NB(@althonos): we use inclusive 1-based ranges in the data model
        # but slicing expects 0-based ranges with exclusive ends
        bgc = self.source[self.start - 1 : self.end]
        bgc.id = bgc.name = self.id
        bgc.seq.alphabet = Bio.Alphabet.generic_dna

        # copy sequence annotations
        bgc.annotations["topology"] = "linear"
        bgc.annotations["organism"] = self.source.annotations.get("organism")
        bgc.annotations["source"] = self.source.annotations.get("source")
        bgc.annotations["comment"] = ["Detected with GECCO v{}".format(__version__)]

        # add proteins as CDS features
        for gene in self.genes:
            # write gene as a /cds GenBank record
            cds = gene.to_seq_feature()
            cds.location += gene.start - self.start
            bgc.features.append(cds)
            # # write domains as /misc_feature annotations
            for domain in gene.protein.domains:
                misc = domain.to_seq_feature(protein_coordinates=False)
                misc.location += cds.location.start
                bgc.features.append(misc)

        # return the complete BGC
        return bgc


@dataclass(frozen=True)
class FeatureTable(Dumpable, Sized):
    """A table storing condensed domain annotations from different genes.
    """

    class Row(NamedTuple):
        sequence_id: str
        protein_id: str
        start: int
        end: int
        strand: str
        domain: str
        hmm: str
        i_evalue: float
        domain_start: int
        domain_end: int
        bgc_probability: Optional[float] = None

    rows: List[Row] = field(default_factory = list)

    @classmethod
    def from_genes(cls, genes: Iterable[Gene]) -> "FeatureTable":
        """Create a new feature table from an iterable of genes.
        """
        rows = []
        for gene in genes:
            for domain in gene.protein.domains:
                rows.append(
                    cls.Row(
                        sequence_id=gene.source.id,
                        protein_id=gene.protein.id,
                        start=gene.start,
                        end=gene.end,
                        strand=gene.strand.sign,
                        domain=domain.name,
                        hmm=domain.hmm,
                        i_evalue=domain.i_evalue,
                        domain_start=domain.start,
                        domain_end=domain.end,
                        bgc_probability=domain.probability
                    )
                )
        return cls(rows)

    def to_genes(self) -> Iterable[Gene]:
        for _, group in itertools.groupby(self, key=operator.attrgetter("protein_id")):
            rows = list(group)
            source = SeqRecord(id=rows[0].sequence_id, seq=Seq("N"*rows[0].end))
            strand = Strand.Coding if rows[0].strand == "+" else Strand.Reverse
            protein = Protein(rows[0].protein_id, seq=None)
            gene = Gene(source, rows[0].start, rows[0].end, strand, protein)
            for row in rows:
                domain = Domain(row.domain, row.domain_start, row.domain_end, row.hmm, row.i_evalue, row.bgc_probability)
                gene.protein.domains.append(domain)
            yield gene

    def __len__(self) -> int:  # noqa: D105
        return len(self.rows)

    def __iter__(self) -> Iterator[Row]:
        return iter(self.rows)

    def dump(self, fh: TextIO, dialect: str = "excel-tab") -> None:
        """Write the feature table in CSV format to the given file.

        Arguments:
            fh (file-like `object`): A writable file-handle opened in text mode
                to write the feature table to.
            dialect (`str`): The CSV dialect to use. See `csv.list_dialects`
                for allowed values.

        """
        writer = csv.writer(fh, dialect=dialect)
        header = list(self.Row.__annotations__)
        writer.writerow(header)
        for row in self:
            writer.writerow([ getattr(row, col) for col in header ])

    @classmethod
    def load(cls, fh: TextIO, dialect: str = "excel-tab") -> "FeatureTable":
        """Load a feature table in CSV format from a file handle in text mode.
        """

        rows = []
        reader = csv.reader(fh, dialect=dialect)
        header = next(reader)

        columns = [
            (header.index(col), col, getattr(ty, "__args__", [ty])[0])
            for col, ty in cls.Row.__annotations__.items()
            if col in header
        ]

        for line in reader:
            raw = { column: ty(line[index]) for index, column, ty in columns }
            rows.append(cls.Row(**raw))
        return cls(rows)


@dataclass(frozen=True)
class ClusterTable(Dumpable, Sized):
    """A table storing condensed information from several clusters.
    """

    class Row(NamedTuple):
        sequence_id: str
        bgc_id: str
        start: int
        end: int
        average_p: float
        max_p: float
        bgc_types: List[str]
        bgc_types_p: List[float]
        proteins: List[str]
        domains: List[str]

    rows: List[Row] = field(default_factory=list)

    @classmethod
    def from_clusters(cls, clusters: Iterable[Cluster]) -> "ClusterTable":
        """Create a new cluster table from an iterable of clusters.
        """
        rows = []
        for cluster in clusters:
            domains = {d.name for g in cluster.genes for d in g.protein.domains}
            rows.append(
                cls.Row(
                    sequence_id=cluster.source.id,
                    bgc_id=cluster.id,
                    start=cluster.start,
                    end=cluster.end,
                    average_p=cluster.average_probability,
                    max_p=cluster.maximum_probability,
                    bgc_types=cluster.types[:],
                    bgc_types_p=cluster.types_probabilities[:],
                    proteins=[gene.protein.id for gene in cluster.genes],
                    domains=sorted(domains),
                )
            )
        return cls(rows)

    def __len__(self) -> int:  # noqa: D105
        return len(self.rows)

    def __iter__(self) -> Iterator[Row]:
        return iter(self.rows)

    def dump(self, fh: TextIO, dialect: str = "excel-tab") -> None:
        """Write the cluster table in CSV format to the given file.

        Arguments:
            fh (file-like `object`): A writable file-handle opened in text mode
                to write the cluster table to.
            dialect (`str`): The CSV dialect to use. See `csv.list_dialects`
                for allowed values.

        """
        writer = csv.writer(fh, dialect=dialect)
        header = list(self.Row.__annotations__)
        writer.writerow(header)
        for row in self:
            line = []
            for col in header:
                value = getattr(row, col)
                if isinstance(value, list):
                    value = ";".join(map(str, value))
                line.append(value)
            writer.writerow(line)

    @classmethod
    def load(cls, fh: TextIO, dialect: str = "excel-tab") -> "ClusterTable":
        reader = csv.reader(fh, dialect=dialect)
        header = next(reader)
        columns = [
            (header.index(col), col, ty)
            for col, ty in cls.Row.__annotations__.items()
        ]

        rows = []
        for line in reader:
            raw = {}
            for index, column, ty in columns:
                if isinstance(ty, type):
                    raw[column] = ty(line[index])
                elif line[index]:
                    ty = ty.__args__[0]
                    raw[column] = list(map(ty, line[index].split(";")))
                else:
                    raw[column] = []
            rows.append(cls.Row(**raw))

        return cls(rows)
