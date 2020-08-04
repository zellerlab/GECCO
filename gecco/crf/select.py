"""Feature selection using Fisher's Exact Test.
"""

import collections
import typing
from typing import Dict, Iterable, Mapping, Optional, Set, Tuple

import fisher
import numpy
from statsmodels.stats.multitest import fdrcorrection

from ..model import Domain, Protein

if typing.TYPE_CHECKING:
    from ..bgc import BGC


def significance_correction(
    significance: Mapping[str, float], method: str,
) -> Dict[str, float]:
    """Perform FDR correction on the ``significance`` dictionary.

    Arguments:
        significance (dict): A dictionary which maps feature names to Fisher
            p-values.
        method (str): The correction method to use. See allowed values in the
            documentation of `statsmodels.stats.multitest.fdrcorrection`.

    Returns:
        dict: A dictionary which maps feature names to corrected p-values.

    Example:
        >>> s = {"A": 0.6, "B": 0.05, "C": 1, "D": 0}
        >>> sorted(significance_correction(s, method="indep").items())
        [('A', 0.7999999999999999), ('B', 0.1), ('C', 1.0), ('D', 0.0)]

    """
    features = sorted(significance, key=significance.__getitem__)
    pvalues = numpy.array([significance[feature] for feature in features])
    _, corrected = fdrcorrection(pvalues, method=method, is_sorted=True)
    return dict(zip(features, corrected))


def fisher_significance(
    data: Iterable[Protein],
    correction_method: Optional[str] = "indep",
    threshold: float = 0.5,
) -> Dict[str, float]:
    r"""Estimate the significance of each domain in the given proteins.

    Example:
        In the following example, we check the significance of three domains
        (*A*, *B* and *C*) on BGC membership for a training set containing
        7 proteins:

        >>> from gecco.model import *
        >>> data = [
        ...     Protein("prot1", _, [
        ...         Domain("A", _, _, _, _, probability=1),
        ...         Domain("B", _, _, _, _, probability=1),
        ...     ]),
        ...     Protein("prot2", _, [
        ...         Domain("A", _, _, _, _, probability=1),
        ...         Domain("B", _, _, _, _, probability=1),
        ...     ]),
        ...     Protein("prot3", _, [
        ...         Domain("A", _, _, _, _, probability=1),
        ...         Domain("B", _, _, _, _, probability=1),
        ...     ]),
        ...     Protein("prot4", _, [Domain("A", _, _, _, _, probability=1)]),
        ...     Protein("prot5", _, [Domain("A", _, _, _, _, probability=1)]),
        ...     Protein("prot6", _, [
        ...         Domain("C", _, _, _, _, probability=0),
        ...         Domain("B", _, _, _, _, probability=0),
        ...     ]),
        ...     Protein("prot7", _, [Domain("C", _, _, _, _, probability=0)]),
        ... ]
        >>> sorted(fisher_significance_new(data).items())
        [('A', 0.071...), ('B', 0.999...), ('C', 0.071...)]

        Since *A* and *C* only appear in BGC and non BGC proteins respectively,
        the p-value for a two-tailed Fisher Exact Test is under 5%, while *B*,
        which appears in half of the BGC proteins and in half of the non-BGC
        proteins, is not significant with regards to the fisher test.

        It's also possible to get the uncorrected values by giving `None`
        instead of a correction method:

        >>> sorted(fisher_significance_new(data, correction_method=None).items())
        [('A', 0.047...), ('B', 0.999...), ('C', 0.047...)]

    """
    # set of all proteins, +proteins grouped by features
    proteins = set(), set()  # type: ignore
    features = collections.defaultdict(set), collections.defaultdict(set)  # type: ignore

    # collect proteins / features for all data tables
    for protein in data:
        for domain in protein.domains:
            in_bgc = domain.probability > threshold
            proteins[in_bgc].add(protein.id)
            features[in_bgc][domain.name].add(protein.id)

    # make the contigency table for each feature
    significance = {}
    for feature in set(features[False]).union(features[True]):
        pvalue = fisher.pvalue(
            len(features[True][feature]),  # with feature, in BGC
            len(proteins[True]) - len(features[True][feature]),  # without feature, in BGC
            len(features[False][feature]),  # with feature, not in BGC
            len(proteins[False]) - len(features[False][feature]),  # without feature, not in BGC
        )
        significance[feature] = pvalue.two_tail

    # perform multiple test correction if needed
    if correction_method is not None:
        significance = significance_correction(significance, correction_method)

    return significance
