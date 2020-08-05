"""Features extraction and sanitization utilities for `~gecco.crf.ClusterCRF`.
"""
import itertools
import typing
from typing import Iterable, List, Sequence, Dict

from ..model import Gene

if typing.TYPE_CHECKING:
    _S = typing.TypeVar("_S", bound=Sequence[Gene])


def extract_features_group(sequence: Iterable[Gene]) -> List[Dict[str, float]]:
    """Extract features at the group/protein level from an iterable of genes.
    """
    # FIXME: currently (v0.3.0 and later) we have to hide proteins missing
    # domain annotations because the model has only been trained with proteins
    # that had domains (since unannotated proteins do not appear in the input
    # feature table used for training)
    #
    #        Fixing this requires retraining the model so that it is aware of
    #  unannotated proteins in the training data, and learns to ignore them.
    #  When it is done, make sure to edit `postprocessing` as well.
    return [{d.name: 1 - d.i_evalue for d in g.protein.domains } for g in sequence if g.protein.domains]

def extract_labels_group(sequence: Iterable[Gene]) -> List[str]:
    """Extract labels at the group/protein level from an iterable of genes.
    """
    return [
        str(int(g.average_probability > 0.5))  # type: ignore
        for g in sequence
        if g.protein.domains
    ]


def annotate_probabilities_group(sequence: "_S", marginals: List[Dict[str, float]]) -> "_S":
    """Annotate genes with marginals obtained from a CRF at the protein level.
    """
    probabilities = iter(marginals)
    for gene in sequence:
        if not gene.protein.domains:
            continue
        p = next(probabilities)["1"]
        gene.protein.domains = [ d.with_probability(p) for d in gene.protein.domains]
    return sequence


def annotate_probabilities_single(sequence: "_S", marginals: List[Dict[str, float]]) -> "_S":
    """Annotate genes with marginals obtained from a CRF at the domain level.
    """
    domains = [ domain for gene in sequence for domain in gene.protein.domains ]
    for domain, probability in itertools.zip_longest(domains, marginals):
        assert domain is not None
        assert probability is not None
        domain.probability = probability["1"]
    return sequence

def extract_labels_single(sequence: Iterable[Gene]) -> List[str]:
    """Extract labels at the domain level from an iterable of genes.
    """
    return [
        str(int(d.probability > 0.5))  # type: ignore
        for g in sequence
        for d in g.protein.domains
    ]


def extract_features_single(sequence: Iterable[Gene]) -> List[Dict[str, float]]:
    """Extract features at the domain level from an iterable of genes.
    """
    return [{d.name: 1 - d.i_evalue} for g in sequence for d in g.protein.domains]
