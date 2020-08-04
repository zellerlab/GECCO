import itertools
from typing import Iterable, List, Dict

from ..model import Gene


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
