"""Data layer classes storing information needed for gene cluster detection.
"""

import collections
import copy
import csv
import datetime
import enum
import functools
import itertools
import math
import operator
import re
import statistics
import typing
from array import array
from collections import OrderedDict
from collections.abc import Sized
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, TextIO, NamedTuple, Union, Iterator, Set

import Bio
import numpy
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, Reference
from Bio.SeqRecord import SeqRecord

from . import __version__
from .interpro import GOTerm
from ._base import Dumpable, Table, _SELF
from ._meta import patch_locale

if typing.TYPE_CHECKING:
    from numpy.typing import NDArray


__all__ = [
    "ClusterType",
    "Strand",
    "Domain",
    "Protein",
    "Gene",
    "Cluster",
    "FeatureTable",
    "ClusterTable"
]


# fmt: off
class ClusterType(object):
    """An immutable storage for the type of a gene cluster.
    """

    def __init__(self, *names: str) -> None:
        """Create a new product type from one or more base types.

        Example:
            >>> t1 = ClusterType()                    # unknown type
            >>> t2 = ClusterType("Polyketide")        # single type
            >>> t3 = ClusterType("Polyketide", "NRP") # multiple types

        """
        self.names = frozenset(names)

    def __repr__(self) -> str:  # noqa: D105
        return "ClusterType({})".format(", ".join(map(repr, sorted(self.names))))

    def __hash__(self) -> int:  # noqa: D105
        return hash(self.names)

    def __eq__(self, other: object) -> bool:  # noqa: D105
        if not isinstance(other, ClusterType):
            return NotImplemented
        return self.names == other.names

    def __bool__(self) -> bool:  # noqa: D105
        return len(self.names) != 0

    def unpack(self) -> List["ClusterType"]:
        """Unpack a composite `ClusterType` into a list of individual types.

        Example:
            >>> ty = ClusterType("Polyketide", "Saccharide")
            >>> ty.unpack()
            [ClusterType('Polyketide'), ClusterType('Saccharide')]

        """
        return [ ClusterType(x) for x in sorted(self.names) ]


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
            is part of a gene cluster, or `None` if no prediction has been 
            made yet.
        cluster_weight (`float`, optional): The weight for this domain, 
            measuring its importance as infered from the training clusters
            by the CRF model.
        go_terms (`list` of `GOTerm`): The Gene Ontology terms
            for this particular domain.
        go_functions (`list` of `GOTerm`): The Gene Ontology term families for 
            this particular domain. Term families are extracted by taking the 
            highest superclasses (excluding the root) of each Gene Ontology 
            term in the ``molecular_function`` namespace associated with this 
            domain.
        qualifiers (`dict`): A dictionary of feature qualifiers that
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
    cluster_weight: Optional[float] = None
    go_terms: List[GOTerm] = field(default_factory=list)
    go_functions: List[GOTerm] = field(default_factory=list)
    qualifiers: Dict[str, List[str]] = field(default_factory=dict)

    def with_probability(self, probability: Optional[float]) -> "Domain":
        """Copy the current domain and assign it a cluster probability.
        """
        return Domain(
            self.name, self.start, self.end, self.hmm, self.i_evalue, self.pvalue,
            probability,
            self.cluster_weight,
            self.go_terms,
            self.go_functions,
            self.qualifiers.copy()
        )

    def with_cluster_weight(self, cluster_weight: Optional[float]) -> "Domain":
        """Copy the current domain and assign it a cluster weight.
        """
        return Domain(
            self.name, self.start, self.end, self.hmm, self.i_evalue, self.pvalue,
            self.probability,
            cluster_weight,
            self.go_terms,
            self.go_functions,
            self.qualifiers.copy()
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
        qualifiers.setdefault("standard_name", [self.name])
        for go_term in self.go_terms:
            qualifiers.setdefault("db_xref", []).append(go_term.accession)
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

    def with_seq(self, seq: Seq) -> "Protein":
        """Copy the current protein and assign it a new sequence.
        """
        return Protein(self.id, seq, copy.deepcopy(self.domains))

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
    qualifiers: Dict[str, List[str]] = field(default_factory=dict)
    _probability: Optional[float] = field(default_factory=lambda: None)

    @property
    def id(self) -> str:
        """`str`: The identifier of the gene (same as the protein identifier).
        """
        return self.protein.id

    @property
    def average_probability(self) -> Optional[float]:
        """`float`: The average of domain probabilities of being in a cluster.
        """
        if self._probability is not None:
            return self._probability
        p = [d.probability for d in self.protein.domains if d.probability is not None]
        return statistics.mean(p) if p else None

    @property
    def maximum_probability(self) -> Optional[float]:
        """`float`: The highest of domain probabilities of being in a cluster.
        """
        if self._probability is not None:
            return self._probability
        p = [d.probability for d in self.protein.domains if d.probability is not None]
        return max(p) if p else None

    # ---

    # NB(@althonos): Color scheme taken from MIBiG.
    #                This is sorted by priority!
    _FUNCTION_PALETTE = {
        # transporter: blue
        "transporter activity": (0x64, 0x95, 0xed),
        "cargo receptor activity": (0x64, 0x95, 0xed),
        "molecular carrier activity": (0x64, 0x95, 0xed),
        # regulatory: green
        "translation regulator activity": (0x2e, 0x8b, 0x56),
        "molecular function regulator activity": (0x2e, 0x8b, 0x56),
        "transcription regulator activity": (0x2e, 0x8b, 0x56),
        "regulation of molecular function": (0x2e, 0x8b, 0x56),
        "general transcription initiation factor activity": (0x2e, 0x8b, 0x56),
        # core biosynthetic: red
        "toxin activity": (0x81, 0x0e, 0x15),
        "catalytic activity": (0x81, 0x0e, 0x15),
        # additional biosynthetic: pink
        "biosynthetic activity": (0xf1, 0x6d, 0x75),
        # non-biosynthetic: olive green
        "non-biosynthetic activity": (0xbd, 0xb7, 0x6b),
        # unknown: grey
        "unknown": (0x80, 0x80, 0x80),
    }

    def to_seq_feature(self, color: bool = True) -> SeqFeature:
        """Convert the gene to a single feature.
        """
        # NB(@althonos): we use inclusive 1-based ranges in the data model
        # but Biopython expects 0-based ranges with exclusive ends
        loc = FeatureLocation(start=self.start, end=self.end+1, strand=int(self.strand))
        qualifiers = dict(self.qualifiers)
        qualifiers.setdefault("locus_tag", [self.protein.id])
        qualifiers.setdefault("translation", [str(self.protein.seq)])

        # NB(@althonos): Attempt to assign a function for the gene based on the
        #                domain content.
        functions = self.functions()
        qualifiers.setdefault("function", sorted(functions))
        if color:
            for k, hexcolor in self._FUNCTION_PALETTE.items():
                if k in functions:
                    break
            else:
                hexcolor = self._FUNCTION_PALETTE["unknown"]
            # EasyFig qualifier
            qualifiers.setdefault("colour", [" ".join(str(x) for x in hexcolor)])
            # SnapGene qualifiers
            qualifiers.setdefault("ApEinfo_fwdcolor", ["#{:02x}{:02x}{:02x}".format(*hexcolor)])
            qualifiers.setdefault("ApEinfo_revcolor", ["#{:02x}{:02x}{:02x}".format(*hexcolor)])

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

    def functions(self) -> Set[str]:
        """Predict the function(s) of the gene from its domain annotations.
        """
        functions = {
            term.name
            for domain in self.protein.domains
            for term in domain.go_functions
        }
        if not functions:
            functions.add("unknown")
        return functions


@dataclass
class Cluster:
    """A sequence of contiguous genes.

    Attributes:
        id (`str`): The identifier of the gene cluster.
        genes (`list` of `~gecco.model.Gene`): A list of the genes belonging
            to this gene cluster.
        types (`gecco.model.ClusterType`): The putative types of this gene 
            cluster, according to similarity in domain composition with 
            curated clusters.
        types_probabilities (`list` of `float`): The probability with which
            each cluster type was identified (same dimension as the ``types``
            attribute).

    """

    id: str
    genes: List[Gene]
    type: ClusterType
    type_probabilities: Dict[str, float]

    def __init__(
        self,
        id: str,
        genes: Optional[List[Gene]] = None,
        type: ClusterType = ClusterType(),
        type_probabilities: Optional[Dict[str, float]] = None,
    ):  # noqa: D107
        self.id = id
        self.genes = genes or list()
        self.type = type
        self.type_probabilities = type_probabilities or dict()

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
        return  statistics.mean(p) if p else None

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
    ) -> "NDArray[numpy.double]":
        """Compute weighted domain composition with respect to ``all_possible``.

        Arguments:
            all_possible (sequence of `str`, optional): A sequence containing
                all domain names to consider when computing domain composition
                for the cluster. If `None` given, then only domains within the
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
            composition of the gene cluster.

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
            return composition / (composition.sum() or 1)
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
        cluster = self.source[self.start - 1 : self.end]
        cluster.id = cluster.name = self.id

        # copy sequence annotations
        cluster.annotations = self.source.annotations.copy()
        cluster.annotations["topology"] = "linear"
        cluster.annotations["molecule_type"] = "DNA"
        with patch_locale("C"):
            cluster.annotations['date'] = now.strftime("%d-%b-%Y").upper()

        biopython_version = tuple(map(int, Bio.__version__.split(".")))
        if biopython_version < (1, 77):
            from Bio import Alphabet

            cluster.seq.alphabet = Alphabet.generic_dna

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
        cluster.annotations.setdefault("references", []).append(ref)

        # add GECCO-specific annotations as a structured comment
        probabilities = {
            f"{key.lower()}_probability":f"{value:.3f}"
            for key, value in self.type_probabilities.items()
        }
        structured_comment = cluster.annotations.setdefault("structured_comment", OrderedDict())
        structured_comment['GECCO-Data'] = {
            "version": f"GECCO v{__version__}",
            "creation_date": now.isoformat(),
            "cluster_type": ";".join(sorted(self.type.names)) or "Unknown",
            **probabilities,
        }

        # add proteins as CDS features
        for gene in self.genes:
            # write gene as a /cds GenBank record
            cds = gene.to_seq_feature()
            cds.location += -self.start
            cluster.features.append(cds)
            # write domains as /misc_feature annotations
            for domain in gene.protein.domains:
                misc = domain.to_seq_feature(protein_coordinates=False)
                if gene.strand == Strand.Coding:
                    misc.location = FeatureLocation(
                        cds.location.start + misc.location.start,
                        cds.location.start + misc.location.end
                    )
                else:
                    misc.location = FeatureLocation(
                        cds.location.end - misc.location.end,
                        cds.location.end - misc.location.start
                    )
                cluster.features.append(misc)

        # return the complete gene cluster record
        return cluster


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
    cluster_probability: List[Optional[float]] = field(default_factory = list)

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
        cluster_probability: Optional[float]

    if typing.TYPE_CHECKING:

        @typing.overload  # type: ignore
        def __getitem__(self, item: int) -> "FeatureTable.Row":  # noqa: D105
            pass

        @typing.overload
        def __getitem__(self, item: slice) -> "FeatureTable.Row":  # noqa: D105
            pass

        def __getitem__(self, item: Union[slice, int]) -> Union["FeatureTable", "FeatureTable.Row"]:   # noqa: D105
            return super().__getitem__(item)  # type: ignore

        def __iter__(self) -> Iterator["FeatureTable.Row"]:  # type: ignore  # noqa: D105
            return super().__iter__()  # type: ignore

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
                table.cluster_probability.append(domain.probability)
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
                domain = Domain(row.domain, row.domain_start, row.domain_end, row.hmm, row.i_evalue, row.pvalue, row.cluster_probability)
                gene.protein.domains.append(domain)
            yield gene

    def __len__(self) -> int:  # noqa: D105
        return len(self.sequence_id)


def _format_product_type(value: "ClusterType") -> str:
    return ";".join(sorted(value.names))


def _parse_product_type(value: str) -> "ClusterType":
    return ClusterType(*map(str.strip, value.split(";"))) if value.strip() else ClusterType()


@dataclass(frozen=True)
class ClusterTable(Table):
    """A table storing condensed information from several clusters.
    """

    sequence_id: List[str] = field(default_factory = list)
    cluster_id: List[str] = field(default_factory = list)
    start: List[int] = field(default_factory = lambda: array("l"))   # type: ignore
    end: List[int] = field(default_factory = lambda: array("l"))     # type: ignore
    average_p: List[Optional[float]] = field(default_factory = list)
    max_p: List[Optional[float]] = field(default_factory = list)

    type: List[ClusterType] = field(default_factory = list)
    type_p: List[Dict[str, float]] = field(default_factory = list)

    proteins: List[Optional[List[str]]] = field(default_factory = list)
    domains: List[Optional[List[str]]] = field(default_factory = list)

    class Row(NamedTuple):
        """A single row in a cluster table.
        """

        sequence_id: str
        cluster_id: str
        start: int
        end: int
        average_p: Optional[float]
        max_p: Optional[float]

        type: ClusterType
        type_p: Dict[ClusterType, float]

        proteins: Optional[List[str]]
        domains: Optional[List[str]]

    _FORMAT_FIELD = {ClusterType: _format_product_type, **Table._FORMAT_FIELD}
    _PARSE_FIELD = {ClusterType: _parse_product_type, **Table._PARSE_FIELD}

    @classmethod
    def from_clusters(cls, clusters: Iterable[Cluster]) -> "ClusterTable":
        """Create a new cluster table from an iterable of clusters.
        """
        table = cls()
        for cluster in clusters:
            table.sequence_id.append(cluster.source.id)
            table.cluster_id.append(cluster.id)
            table.start.append(cluster.start)
            table.end.append(cluster.end)
            table.average_p.append(cluster.average_probability)
            table.max_p.append(cluster.maximum_probability)

            table.type.append(cluster.type)
            table.type_p.append(cluster.type_probabilities.copy())

            table.proteins.append([ gene.protein.id for gene in cluster.genes ])
            domains = {d.name for g in cluster.genes for d in g.protein.domains}
            table.domains.append(sorted(domains))
        return table

    def __len__(self) -> int:  # noqa: D105
        return len(self.sequence_id)

    if typing.TYPE_CHECKING:

        @typing.overload  # type: ignore
        def __getitem__(self, item: int) -> "ClusterTable.Row":  # noqa: D105
            pass

        @typing.overload
        def __getitem__(self, item: slice) -> "ClusterTable.Row":  # noqa: D105
            pass

        def __getitem__(self, item: Union[slice, int]) -> Union["ClusterTable", "ClusterTable.Row"]:   # noqa: D105
            return super().__getitem__(item)  # type: ignore

        def __iter__(self) -> Iterator["ClusterTable.Row"]:  # type: ignore  # noqa: D105
            return super().__iter__()  # type: ignore

    def dump(self, fh: TextIO, dialect: str = "excel-tab", header: bool = True) -> None:  # noqa: D102
        writer = csv.writer(fh, dialect=dialect)
        column_names = list(self.__annotations__)
        type_p_col = column_names.index("type_p")
        column_names.pop(type_p_col)
        optional = self._optional_columns()

        # do not write optional columns if they are completely empty
        for name in optional:
            if all(x is None for x in getattr(self, name)):
                column_names.remove(name)

        # build column formatters
        columns = [getattr(self, name) for name in column_names]
        formatters = [self._FORMAT_FIELD[self.Row.__annotations__[name]] for name in column_names]

        # add probability columns if any row countains probabilities
        if any(d for d in self.type_p):
            classes = sorted({
                name
                for probabilities in self.type_p
                for name in probabilities.keys()
            })
            for name in classes:
                column_names.insert(type_p_col, f"{name.lower()}_probability")
                columns.insert(type_p_col, [self.type_p[i][name] for i in range(len(self))])
                formatters.insert(type_p_col, str)  # type: ignore
                type_p_col += 1

         # write header if desired
        if header:
            writer.writerow(column_names)
        # write each row
        for i in range(len(self)):
            writer.writerow([format(col[i]) for col,format in zip(columns, formatters)])

    @classmethod
    def load(cls, fh: TextIO, dialect: str = "excel-tab") -> "ClusterTable":  # noqa: D102
        table = cls()
        reader = csv.reader(fh, dialect=dialect)
        header = next(reader)

        # check that if a column is missing, it is one of the optional values
        optional = cls._optional_columns()
        missing = set(cls.__annotations__).difference({"type_p"}).difference(header)
        missing_required = missing.difference(optional)
        if missing_required:
            raise ValueError("table is missing columns: {}".format(", ".join(missing_required)))

        # get the name of each column and check which columns are optional
        columns = [getattr(table, col, None) for col in header]
        parsers = [
            cls._PARSE_FIELD[table.Row.__annotations__[col]]
            if col in table.Row.__annotations__
            else None
            for col in header
        ]

        # extract elements from the CSV rows
        for row in reader:
            probabilities = {}
            for name, col, value, parser in zip(header, columns, row, parsers):
                if parser is not None and col is not None:
                    col.append(parser(value))
                elif name.endswith("_probability"):
                    probabilities[name[:-12]] = float(value)
            table.type_p.append(probabilities)

        for col in missing:
            getattr(table, col).extend(None for _ in range(len(table)))
        return table


@dataclass(frozen=True)
class GeneTable(Table):
    """A table storing gene coordinates and optional cluster probabilities.
    """

    sequence_id: List[str] = field(default_factory = list)
    protein_id: List[str] = field(default_factory = list)
    start: List[int] = field(default_factory = lambda: array("l"))    # type: ignore
    end: List[int] = field(default_factory = lambda: array("l"))      # type: ignore
    strand: List[str] = field(default_factory = list)
    average_p: List[Optional[float]] = field(default_factory = list)
    max_p: List[Optional[float]] = field(default_factory = list)

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

    if typing.TYPE_CHECKING:

        @typing.overload  # type: ignore
        def __getitem__(self, item: int) -> "GeneTable.Row":  # noqa: D105
            pass

        @typing.overload
        def __getitem__(self, item: slice) -> "GeneTable.Row":  # noqa: D105
            pass

        def __getitem__(self, item: Union[slice, int]) -> Union["GeneTable", "GeneTable.Row"]:   # noqa: D105
            return super().__getitem__(item)  # type: ignore

        def __iter__(self) -> Iterator["GeneTable.Row"]:  # type: ignore  # noqa: D105
            return super().__iter__()  # type: ignore

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
        for row in typing.cast(Iterable["GeneTable.Row"], self):
            source = SeqRecord(id=row.sequence_id, seq=_UnknownSeq())
            strand = Strand.Coding if row.strand == "+" else Strand.Reverse
            seq = Seq("X" * (row.end - row.start // 3))
            protein = Protein(row.protein_id, seq=_UnknownSeq())
            yield Gene(source, row.start, row.end, strand, protein, _probability=row.average_p)

    def __len__(self) -> int:  # noqa: D105
        return len(self.protein_id)
