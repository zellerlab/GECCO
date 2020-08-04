"""Features extraction and sanitization utilities for `~gecco.crf.ClusterCRF`.
"""
import itertools
from typing import Iterable, List, Dict

from ..model import Gene


def extract_features_group(sequence: Iterable[Gene]) -> List[Dict[str, float]]:
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
    return [str(int(g.average_probability > 0.5)) for g in sequence if g.protein.domains]


def annotate_probabilities_group(
    sequence: Iterable[Gene],
    marginals: List[List[Dict[str, float]]]
) -> List[Dict[str, float]]:
    probabilities = iter(marginals)
    for gene in sequence:
        if not gene.protein.domains:
            continue
        p = next(probabilities)["1"]
        gene.protein.domains = [ d.with_probability(p) for d in gene.protein.domains]
    return sequence


def annotate_probabilities_single(
    sequence: Iterable[Gene],
    marginals: List[List[Dict[str, float]]]
) -> List[Dict[str, float]]:
    domains = [ domain for gene in sequence for domain in gene.protein.domains ]
    for domain, probability in itertools.zip_longest(domains, marginals):
        assert domain is not None
        assert probability is not None
        domain.probability = probability["1"]
    return sequence

def extract_labels_single(sequence: Iterable[Gene]) -> List[str]:
    return [str(int(d.probability > 0.5)) for g in sequence for d in g.protein.domains]


def extract_features_single(sequence: Iterable[Gene]) -> List[Dict[str, float]]:
    return [{d.name: 1 - d.i_evalue} for g in sequence for d in g.protein.domains]
