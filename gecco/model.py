"""Data layer classes storing information needed for BGC detection.
"""

import collections
import csv
import datetime
import enum
import functools
import itertools
import math
import operator
import re
import typing
from array import array
from collections import OrderedDict
from collections.abc import Sized
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, TextIO, NamedTuple, Union, Iterator

import Bio
import numpy
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, Reference
from Bio.SeqRecord import SeqRecord

from . import __version__
from ._base import Dumpable, Table
from ._meta import patch_locale


__all__ = [
    "ProductType",
    "Strand",
    "Domain",
    "Protein",
    "Gene",
    "Cluster",
    "FeatureTable",
    "ClusterTable"
]


# fmt: off
class ProductType(enum.IntFlag):
    """A flag to declare the type of product synthesized by a gene cluster.
    """

    Unknown    = 0b00000000
    Alkaloid   = 0b00000010
    Polyketide = 0b00000100
    RiPP       = 0b00001000
    Saccharide = 0b00010000
    Terpene    = 0b00100000
    NRP        = 0b01000000

    @classmethod
    def pack(cls, members: Iterable["ProductType"]) -> "ProductType":
        """Pack together a list of individual product types.

        Example:
            >>> types = [ProductType.Polyketide, ProductType.Saccharide]
            >>> ProductType.pack(types)
            <ProductType.Saccharide|Polyketide: 20>

        """
        return functools.reduce(operator.or_, members, cls.Unknown)

    def unpack(self) -> List["ProductType"]:
        """Unpack a composite `ProductType` into a list of individual types.

        Example:
            >>> ty = ProductType.Polyketide | ProductType.Saccharide
            >>> ty.unpack()
            [<ProductType.Polyketide: 4>, <ProductType.Saccharide: 16>]

        """
        return [ x for x in ProductType.__members__.values() if (x & self) ]


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
        pvalue (`float`): The p-value reported by ``hmmsearch`` that measure
            how likely the domain score is.
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
    pvalue: float
    probability: Optional[float] = None
    qualifiers: Mapping[str, Union[str, List[str]]] = field(default_factory=dict)

    def with_probability(self, probability: Optional[float]) -> "Domain":
        """Copy the current domain and assign it a BGC probability.
        """
        return Domain(
            self.name, self.start, self.end, self.hmm, self.i_evalue, self.pvalue,
            probability, self.qualifiers.copy()
        )

    def to_seq_feature(self, protein_coordinates: bool = False) -> SeqFeature:
        """Convert the domain to a single feature.

        Arguments:
            protein_coordinates (`bool`): Set to `True` for the feature
                coordinates to be given in amino-acids, or to `False` in
                nucleotides.

        """
        stride = 1 if protein_coordinates else 3
        loc = FeatureLocation((self.start-1)*stride, self.end*stride)
        qualifiers = dict(self.qualifiers)
        qualifiers.setdefault("standard_name", self.name)
        return SeqFeature(location=loc, type="misc_feature", qualifiers=qualifiers)


@dataclass(frozen=True)
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
        record = SeqRecord(self.seq, id=self.id, name=self.id)

        biopython_version = tuple(map(int, Bio.__version__.split(".")))
        if biopython_version < (1, 77):
            from Bio import Alphabet

            record.seq.alphabet = Alphabet.generic_protein

        return record

    def with_domains(self, domains: Iterable[Domain]) -> "Protein":
        """Copy the current protein and assign it new domains.
        """
        return Protein(self.id, self.seq, list(domains))


@dataclass(frozen=True)
class Gene:
    """A nucleotide sequence coding a protein.

    Attributes:
        source (`~Bio.SeqRecord.SeqRecord`): The DNA sequence this gene was
            found in, as a Biopython record.
        start (`int`): The index of the leftmost nucleotide of the gene within
            the source sequence, independent of the strandedness.
        end (`int`): The index of the rightmost nucleotide of the gene within
            the source sequence.
        strand (`~gecco.model.Strand`): The strand where the gene is located.
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
    _probability: Optional[float] = field(default_factory=lambda: None)

    @property
    def id(self) -> str:
        """`str`: The identifier of the gene (same as the protein identifier).
        """
        return self.protein.id

    @property
    def average_probability(self) -> Optional[float]:
        """`float`: The average of domain probabilities of being biosynthetic.
        """
        if self._probability is not None:
            return self._probability
        p = [d.probability for d in self.protein.domains if d.probability is not None]
        return sum(p) / len(p) if p else None

    @average_probability.setter
    def average_probability(self, probability: Optional[float]):
        self._probability = probability

    @property
    def maximum_probability(self) -> Optional[float]:
        """`float`: The highest of domain probabilities of being biosynthetic.
        """
        if self._probability is not None:
            return self._probability
        p = [d.probability for d in self.protein.domains if d.probability is not None]
        return max(p) if p else None

    # ---

    def to_seq_feature(self) -> SeqFeature:
        """Convert the gene to a single feature.
        """
        # NB(@althonos): we use inclusive 1-based ranges in the data model
        # but Biopython expects 0-based ranges with exclusive ends
        loc = FeatureLocation(start=self.start, end=self.end+1, strand=int(self.strand))
        qualifiers = dict(self.qualifiers)
        qualifiers.setdefault("locus_tag", self.protein.id)
        qualifiers.setdefault("translation", str(self.protein.seq))
        return SeqFeature(location=loc, type="CDS", qualifiers=qualifiers)

    def with_protein(self, protein: "Protein") -> "Gene":
        """Copy the current gene and assign it a different protein.
        """
        return Gene(
            self.source, self.start, self.end, self.strand, protein,
            self.qualifiers.copy(),
            _probability=self._probability,
        )

    def with_source(self, source: "SeqRecord") -> "Gene":
        """Copy the current gene and assign it a different source.
        """
        return Gene(
            source, self.start, self.end, self.strand, self.protein,
            self.qualifiers.copy(),
            _probability=self._probability,
        )

    def with_probability(self, probability: float) -> "Gene":
        """Copy the current gene and assign it a different probability.
        """
        return Gene(
            self.source, self.start, self.end, self.strand,
            self.protein.with_domains([
                domain.with_probability(probability)
                for domain in self.protein.domains
            ]),
            self.qualifiers.copy(),
            _probability=probability
        )


@dataclass
class Cluster:
    """A sequence of contiguous genes with biosynthetic activity.

    Attributes:
        id (`str`): The identifier of the gene cluster.
        genes (`list` of `~gecco.model.Gene`): A list of the genes belonging
            to this gene cluster.
        types (`gecco.model.ProductType`): The putative types of product
            synthesized by this gene cluster, according to similarity in
            domain composition with curated clusters.
        types_probabilities (`list` of `float`): The probability with which
            each BGC type was identified (same dimension as the ``types``
            attribute).

    """

    id: str
    genes: List[Gene]
    type: ProductType
    type_probabilities: Mapping[ProductType, float]

    def __init__(
        self,
        id: str,
        genes: Optional[List[Gene]] = None,
        type: ProductType = ProductType.Unknown,
        type_probabilities: Optional[Dict[ProductType, float]] = None,
    ):  # noqa: D107
        self.id = id
        self.genes = genes or list()
        self.type = type
        self.type_probabilities = type_probabilities or dict()

        # if len(self.types) != len(self.types_probabilities):
        #     err = "type and type probability lists must have the same dimensions"
        #     raise ValueError(err)


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
        self,
        all_possible: Optional[Sequence[str]] = None,
        normalize: bool = True,
        minlog_weights: bool = False,
        pvalue: bool = True,
    ) -> numpy.ndarray:
        """Compute weighted domain composition with respect to ``all_possible``.

        Arguments:
            all_possible (sequence of `str`, optional): A sequence containing
                all domain names to consider when computing domain composition
                for the BGC. If `None` given, then only domains within the
                cluster are taken into account.
            normalize (`bool`): Normalize the composition vector so that it
                sums to 1.
            minlog_weights (`bool`): Compute weight for each domain as
                :math:`-log_10(v)` (where :math:`v` is either the ``pvalue``
                or the ``i_evalue``, depending on the value of ``normalize``).
                Use :math:`1 - v` otherwise.
            pvalue (`bool`): Compute composition weights using the ``pvalue``
                of each domain, instead of the ``i_evalue``.

        Returns:
            `~numpy.ndarray`: A numerical array containing the relative domain
            composition of the BGC.

        """
        domains = [d for gene in self.genes for d in gene.protein.domains]
        names = numpy.array([domain.name for domain in domains])
        field = operator.attrgetter("pvalue" if pvalue else "i_evalue")
        if minlog_weights:
            weights = numpy.array([- math.log10(field(domain)) for domain in domains])
        else:
            weights = numpy.array([1 - field(domain) for domain in domains])

        unique_names = set(names)
        if all_possible is None:
            all_possible = numpy.unique(names)
        composition = numpy.zeros(len(all_possible))
        for i, dom in enumerate(all_possible):
            if dom in unique_names:
                composition[i] = numpy.sum(weights[names == dom])
        if normalize:
            return composition / (composition.sum() or 1)  # type: ignore
        return composition

    # ---

    def to_seq_record(self) -> SeqRecord:
        """Convert the cluster to a single record.

        Annotations of the source sequence are kept intact if they don't
        overlap with the cluster boundaries. Component genes are added on the
        record as *CDS* features. Annotated protein domains are added as
        *misc_feature*.

        """
        # store time of record creation
        now = datetime.datetime.now()

        # NB(@althonos): we use inclusive 1-based ranges in the data model
        # but slicing expects 0-based ranges with exclusive ends
        bgc = self.source[self.start - 1 : self.end]
        bgc.id = bgc.name = self.id

        # copy sequence annotations
        bgc.annotations = self.source.annotations.copy()
        bgc.annotations["topology"] = "linear"
        bgc.annotations["molecule_type"] = "DNA"
        with patch_locale("C"):
            bgc.annotations['date'] = now.strftime("%d-%b-%Y").upper()

        biopython_version = tuple(map(int, Bio.__version__.split(".")))
        if biopython_version < (1, 77):
            from Bio import Alphabet

            bgc.seq.alphabet = Alphabet.generic_dna

        # add GECCO preprint as a reference
        ref = Reference()
        ref.title = "Accurate de novo identification of biosynthetic gene clusters with GECCO"
        ref.journal = "bioRxiv (2021.05.03.442509)"
        ref.comment = "doi:10.1101/2021.05.03.442509"
        ref.authors = ", ".join([
            "Laura M Carroll",
            "Martin Larralde",
            "Jonas Simon Fleck",
            "Ruby Ponnudurai",
            "Alessio Milanese",
            "Elisa Cappio Barazzone",
            "Georg Zeller"
        ])
        bgc.annotations.setdefault("references", []).append(ref)

        # add GECCO-specific annotations as a structured comment
        structured_comment = bgc.annotations.setdefault("structured_comment", OrderedDict())
        structured_comment['GECCO-Data'] = {
            "version": f"GECCO v{__version__}",
            "creation_date": now.isoformat(),
            "biosyn_class": ";".join(sorted(ty.name for ty in self.type.unpack())) or "Unknown",
            "alkaloid_probability": self.type_probabilities.get(ProductType.Alkaloid, 0.0),
            "polyketide_probability": self.type_probabilities.get(ProductType.Polyketide, 0.0),
            "ripp_probability": self.type_probabilities.get(ProductType.RiPP, 0.0),
            "saccharide_probability": self.type_probabilities.get(ProductType.Saccharide, 0.0),
            "terpene_probability": self.type_probabilities.get(ProductType.Terpene, 0.0),
            "nrp_probability": self.type_probabilities.get(ProductType.NRP, 0.0),
        }

        # add proteins as CDS features
        for gene in self.genes:
            # write gene as a /cds GenBank record
            cds = gene.to_seq_feature()
            cds.location += -self.start
            bgc.features.append(cds)
            # write domains as /misc_feature annotations
            for domain in gene.protein.domains:
                misc = domain.to_seq_feature(protein_coordinates=False)
                misc.location += cds.location.start
                bgc.features.append(misc)

        # return the complete BGC
        return bgc


class _UnknownSeq(Seq):
    """An unknown sequence that uses limited memory.

    Used by `FeatureTable.to_genes` to fake the `Gene.source.seq` attribute.
    """

    def __init__(self) -> None:
        super().__init__(data="")

    @typing.overload
    def __getitem__(self, index: int) -> str:
        pass

    @typing.overload
    def __getitem__(self, index: slice) -> Seq:
        pass

    def __getitem__(self, index: Union[slice, int]) -> Union[str, Seq]:
        if isinstance(index, slice):
            return Seq("N" * ((index.stop - index.start) // (index.step or 1)) )
        return "N"


@dataclass(frozen=True)
class FeatureTable(Table):
    """A table storing condensed domain annotations from different genes.
    """

    sequence_id: List[str] = field(default_factory = list)
    protein_id: List[str] = field(default_factory = list)
    start: List[int] = field(default_factory = lambda: array("l"))          # type: ignore
    end: List[int] = field(default_factory = lambda: array("l"))            # type: ignore
    strand: List[str] = field(default_factory = list)
    domain: List[str] = field(default_factory = list)
    hmm: List[str] = field(default_factory = list)
    i_evalue: List[float] = field(default_factory = lambda: array("d"))     # type: ignore
    pvalue: List[float] = field(default_factory = lambda: array("d"))       # type: ignore
    domain_start: List[int] = field(default_factory = lambda: array("l"))   # type: ignore
    domain_end: List[int] = field(default_factory = lambda: array("l"))     # type: ignore
    bgc_probability: List[Optional[float]] = field(default_factory = list)

    class Row(NamedTuple):
        """A single row in a feature table.
        """

        sequence_id: str
        protein_id: str
        start: int
        end: int
        strand: str
        domain: str
        hmm: str
        i_evalue: float
        pvalue: float
        domain_start: int
        domain_end: int
        bgc_probability: Optional[float]

    @classmethod
    def from_genes(cls, genes: Iterable[Gene]) -> "FeatureTable":
        """Create a new feature table from an iterable of genes.
        """
        table = cls()
        for gene in genes:
            for domain in gene.protein.domains:
                table.sequence_id.append(gene.source.id)
                table.protein_id.append(gene.protein.id)
                table.start.append(gene.start)
                table.end.append(gene.end)
                table.strand.append(gene.strand.sign)
                table.domain.append(domain.name)
                table.hmm.append(domain.hmm)
                table.i_evalue.append(domain.i_evalue)
                table.pvalue.append(domain.pvalue)
                table.domain_start.append(domain.start)
                table.domain_end.append(domain.end)
                table.bgc_probability.append(domain.probability)
        return table

    def to_genes(self) -> Iterable[Gene]:
        """Convert a feature table to actual genes.

        Since the source sequence cannot be known, a *dummy* sequence is
        built for each gene of size ``gene.end``, so that each gene can still
        be converted to a `~Bio.SeqRecord.SeqRecord` if needed.
        """
        # group rows by protein/gene ID
        protein_indices = collections.defaultdict(list)
        for i, protein_id in enumerate(self.protein_id):
            protein_indices[protein_id].append(i)
        # yield genes in order
        for protein_id in sorted(protein_indices):
            rows = [self[i] for i in protein_indices[protein_id]]
            assert all(x.sequence_id == rows[0].sequence_id for x in rows)
            assert all(x.protein_id == rows[0].protein_id for x in rows)
            assert all(x.start == rows[0].start for x in rows)
            assert all(x.end == rows[0].end for x in rows)
            source = SeqRecord(id=rows[0].sequence_id, seq=_UnknownSeq())
            strand = Strand.Coding if rows[0].strand == "+" else Strand.Reverse
            protein = Protein(rows[0].protein_id, seq=_UnknownSeq)
            gene = Gene(source, rows[0].start, rows[0].end, strand, protein)
            for row in rows:
                domain = Domain(row.domain, row.domain_start, row.domain_end, row.hmm, row.i_evalue, row.pvalue, row.bgc_probability)
                gene.protein.domains.append(domain)
            yield gene

    def __len__(self) -> int:  # noqa: D105
        return len(self.sequence_id)


def _format_product_type(value: "ProductType") -> str:
    types = value.unpack() or [ProductType.Unknown]
    return ";".join(sorted(map(operator.attrgetter("name"), types)))


def _parse_product_type(value: str) -> "ProductType":
    types = [ProductType.__members__[x] for x in value.split(";")]
    return ProductType.pack(types)


@dataclass(frozen=True)
class ClusterTable(Table):
    """A table storing condensed information from several clusters.
    """

    sequence_id: List[str] = field(default_factory = list)
    bgc_id: List[str] = field(default_factory = list)
    start: List[int] = field(default_factory = lambda: array("l"))          # type: ignore
    end: List[int] = field(default_factory = lambda: array("l"))            # type: ignore
    average_p: List[Optional[float]] = field(default_factory = list)    # type: ignore
    max_p: List[Optional[float]] = field(default_factory = list)        # type: ignore

    type: List[ProductType] = field(default_factory = list)
    alkaloid_probability: List[Optional[float]] = field(default_factory = list)    # type: ignore
    polyketide_probability: List[Optional[float]] = field(default_factory = list)  # type: ignore
    ripp_probability: List[Optional[float]] = field(default_factory = list)        # type: ignore
    saccharide_probability: List[Optional[float]] = field(default_factory = list)  # type: ignore
    terpene_probability: List[Optional[float]] = field(default_factory = list)     # type: ignore
    nrp_probability: List[Optional[float]] = field(default_factory = list)         # type: ignore

    proteins: List[Optional[List[str]]] = field(default_factory = list)
    domains: List[Optional[List[str]]] = field(default_factory = list)

    class Row(NamedTuple):
        """A single row in a cluster table.
        """

        sequence_id: str
        bgc_id: str
        start: int
        end: int
        average_p: Optional[float]
        max_p: Optional[float]
        type: ProductType
        alkaloid_probability: Optional[float]
        polyketide_probability: Optional[float]
        ripp_probability: Optional[float]
        saccharide_probability: Optional[float]
        terpene_probability: Optional[float]
        nrp_probability: Optional[float]
        proteins: Optional[List[str]]
        domains: Optional[List[str]]

    _FORMAT_FIELD = {ProductType: _format_product_type, **Table._FORMAT_FIELD}
    _PARSE_FIELD = {ProductType: _parse_product_type, **Table._PARSE_FIELD}

    @classmethod
    def from_clusters(cls, clusters: Iterable[Cluster]) -> "ClusterTable":
        """Create a new cluster table from an iterable of clusters.
        """
        table = cls()
        for cluster in clusters:
            table.sequence_id.append(cluster.source.id)
            table.bgc_id.append(cluster.id)
            table.start.append(cluster.start)
            table.end.append(cluster.end)
            table.average_p.append(cluster.average_probability)
            table.max_p.append(cluster.maximum_probability)

            table.type.append(cluster.type)
            table.alkaloid_probability.append(cluster.type_probabilities.get(ProductType.Alkaloid, 0))
            table.polyketide_probability.append(cluster.type_probabilities.get(ProductType.Polyketide, 0))
            table.ripp_probability.append(cluster.type_probabilities.get(ProductType.RiPP, 0))
            table.saccharide_probability.append(cluster.type_probabilities.get(ProductType.Saccharide, 0))
            table.terpene_probability.append(cluster.type_probabilities.get(ProductType.Terpene, 0))
            table.nrp_probability.append(cluster.type_probabilities.get(ProductType.NRP, 0))

            table.proteins.append([ gene.protein.id for gene in cluster.genes ])
            domains = {d.name for g in cluster.genes for d in g.protein.domains}
            table.domains.append(sorted(domains))
        return table

    def __len__(self) -> int:  # noqa: D105
        return len(self.sequence_id)


@dataclass(frozen=True)
class GeneTable(Table):
    """A table storing gene coordinates and optional biosynthetic probabilities.
    """

    sequence_id: List[str] = field(default_factory = list)
    protein_id: List[str] = field(default_factory = list)
    start: List[int] = field(default_factory = lambda: array("l"))    # type: ignore
    end: List[int] = field(default_factory = lambda: array("l"))      # type: ignore
    strand: List[str] = field(default_factory = list)
    average_p: List[Optional[float]] = field(default_factory = list)  # type: ignore
    max_p: List[Optional[float]] = field(default_factory = list)      # type: ignore

    class Row(NamedTuple):
        """A single row in a gene table.
        """

        sequence_id: str
        protein_id: str
        start: int
        end: int
        strand: str
        average_p: Optional[float]
        max_p: Optional[float]

    @classmethod
    def from_genes(cls, genes: Iterable[Gene]) -> "GeneTable":
        """Create a new gene table from an iterable of genes.
        """
        table = cls()
        for gene in genes:
            table.sequence_id.append(gene.source.id)
            table.protein_id.append(gene.protein.id)
            table.start.append(gene.start)
            table.end.append(gene.end)
            table.strand.append(gene.strand.sign)
            table.average_p.append(gene.average_probability)
            table.max_p.append(gene.maximum_probability)
        return table

    def to_genes(self) -> Iterable[Gene]:
        """Convert a gene table to actual genes.

        Since the source sequence cannot be known, a *dummy* sequence is
        built for each gene of size ``gene.end``, so that each gene can still
        be converted to a `~Bio.SeqRecord.SeqRecord` if needed.

        """
        for row in self:
            source = SeqRecord(id=row.sequence_id, seq=_UnknownSeq())
            strand = Strand.Coding if row.strand == "+" else Strand.Reverse
            seq = Seq("X" * (row.end - row.start // 3))
            protein = Protein(row.protein_id, seq=_UnknownSeq())
            yield Gene(source, row.start, row.end, strand, protein, _probability=row.average_p)

    def __len__(self) -> int:  # noqa: D105
        return len(self.protein_id)
