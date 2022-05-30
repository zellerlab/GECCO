"""Features extraction and sanitization utilities for `~gecco.crf.ClusterCRF`.
"""
import itertools
import typing
from typing import Iterable, List, Sequence, Dict

from ..model import Gene

if typing.TYPE_CHECKING:
    _S = typing.TypeVar("_S", bound=Sequence[Gene])


def extract_features_protein(
    sequence: Iterable[Gene], empty: bool = True
) -> List[Dict[str, bool]]:
    """Extract features at the gene level.

    If several domains are part of the same gene, they are grouped at the
    same position. Genes without domain annotation are represented with an
    empty dictionary.

    Arguments:
        sequence (iterable of `~gecco.model.Gene`): The genes to extract the
            features from.
        empty (`bool`): Whether or not to skip empty genes. With `False`,
            no empty feature dictionary will be emitted, and genes will be
            skipped. With `True`, an empty feature dictionary will be emitted
            for every gene without domain annotation.

    """
    return [
        {domain.name: True for domain in gene.protein.domains}
        for gene in sequence
        if gene.protein.domains or empty
    ]


def extract_features_domain(
    sequence: Iterable[Gene], empty: bool = True
) -> List[Dict[str, bool]]:
    """Extract features at the domain level."""
    features: List[Dict[str, bool]] = []
    for gene in sequence:
        if gene.protein.domains:
            features.extend({domain.name: True} for domain in gene.protein.domains)
        elif empty:
            features.append({})
    return features


def extract_labels_protein(sequence: Iterable[Gene], empty: bool = True) -> List["str"]:
    """Extract labels at the gene level for training."""
    return [
        "1" if typing.cast(float, gene.average_probability) > 0.5 else "0"
        for gene in sequence
        if gene.protein.domains or empty
    ]


def extract_labels_domain(sequence: Iterable[Gene], empty: bool = True) -> List["str"]:
    """Extract labels at the domain level for training."""
    labels: List[str] = []
    for gene in sequence:
        if gene.protein.domains:
            labels.extend(
                "1" if typing.cast(float, domain.probability) > 0.5 else "0"
                for domain in gene.protein.domains
            )
        elif empty:
            labels.append("1" if typing.cast(float, gene.average_probability) > 0.5 else "0")
    return labels


def annotate_probabilities_protein(
    sequence: Sequence[Gene],
    probabilities: Sequence[float],
    empty: bool = True,
) -> Iterable[Gene]:
    """Annotate genes with marginals obtained from a CRF at the protein level.

    Arguments:
        sequence (iterable of `~gecco.model.Gene`): The genes to annotate the
            features from.
        probabilities (iterable of `float`): The biosynthetic probabilities
            that have been computed by the CRF.
        empty (`bool`): Whether or not to empty genes where skipped. With
            `False`, no empty feature dictionary have been emitted, so empty
            genes have no probability. With `True`, empty genes have a
            probability that can be extracted.

    """
    genes = [gene for gene in sequence if gene.protein.domains or empty]
    if len(genes) != len(probabilities):
        raise ValueError("gene and probability lists don't have the same length")
    for gene, probability in zip(genes, probabilities):
        yield gene.with_probability(probability)


def annotate_probabilities_domain(
    sequence: Iterable[Gene],
    probabilities: Iterable[float],
    empty: bool = True,
) -> Iterable[Gene]:
    """Annotate genes with marginals obtained from a CRF at the domain level."""
    genes = iter(sequence)
    probas = iter(probabilities)
    for gene in genes:
        if gene.protein.domains:
            yield gene.with_protein(
                gene.protein.with_domains(
                    [
                        domain.with_probability(p)
                        for domain, p in zip(gene.protein.domains, probas)
                    ]
                )
            )
        elif empty:
            yield gene.with_probability(next(probas))
    if next(probas, None) is not None:
        raise ValueError("gene and probability lists don't have the same length")
