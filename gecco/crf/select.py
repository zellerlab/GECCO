"""Feature selection using Fisher's Exact Test.
"""

import collections
import typing
from typing import Dict, Iterable, Mapping, Optional, Set, Tuple

import numpy

from ..model import Domain, Protein
from .._meta import requires

if typing.TYPE_CHECKING:
    from ..bgc import BGC

_CORRECTION_METHODS = {
    "bonferroni",
    "sidak",
    "holm-sidak",
    "holm",
    "simes-hochberg",
    "hommel",
    "fdr_bh",
    "fdr_by",
    "fdr_tsbh",
    "fdr_tsbky",
}

@requires("statsmodels.stats.multitest")
def significance_correction(
    significance: Mapping[str, float], method: str,
) -> Dict[str, float]:
    """Perform FDR correction on the ``significance`` dictionary.

    Arguments:
        significance (dict): A dictionary which maps feature names to Fisher
            p-values.
        method (str): The correction method to use. See allowed values in the
            documentation of `statsmodels.stats.multitest.multipletests`.

    Returns:
        dict: A dictionary which maps feature names to corrected p-values.

    Example:
        >>> s = {"A": 0.6, "B": 0.05, "C": 1, "D": 0}
        >>> sorted(significance_correction(s, method="fdr_bh").items())
        [('A', 0.7999999999999999), ('B', 0.1), ('C', 1.0), ('D', 0.0)]

    """
    features = sorted(significance, key=significance.__getitem__)
    pvalues = numpy.array([significance[feature] for feature in features])
    _, corrected, _, _ = multitest.multipletests(pvalues, method=method, is_sorted=True)
    return dict(zip(features, corrected))


@requires("fisher")
def fisher_significance(
    proteins: Iterable[Protein],
    correction_method: Optional[str] = "fdr_bh",
) -> Dict[str, float]:
    r"""Estimate the significance of each domain in the given proteins.

    For each feature $F$, we create the following contingency table, by
    counting how many time the feature appears or does not appear in BGC
    proteins, and non-BGC proteins:

    +-------------+-------------------------+-------------------------------+
    | *proteins*  |    with feature $F$     |      without feature $F $     |
    +=============+=========================+===============================+
    | within BGC  |    :math:`N_{F,BGC}`    |    :math:`N_{\bar{F},BGC}`    |
    +-------------+-------------------------+-------------------------------+
    | outside BGC | :math:`N_{F,\bar{BGC}}` | :math:`N_{\bar{F},\bar{BGC}}` |
    +-------------+-------------------------+-------------------------------+

    Then, we run a Fisher Exact Test on this distribution, which gives us the
    probability to observe the same table under the hypothesis of independece
    of the two variables.

    Arguments:
        proteins (iterable of `~gecco.model.Protein`): An iterable yielding
            annotated proteins which domains to estimate the significance of.
            **Domains must have a ``probability`` of 1 if they are part of
            BGC, or of 0 if they are not.**
        correction_method (`str`, optional): The name of the multiple test
            correction method to use when computing significance, or `None` to
            skip correction. See `statsmodels.stats.multitest.multipletests`
            for allowed values.

    Returns:
        `dict`: A dictionary which to each feature associates the p-value of
        the two-tail Fisher Exact Test in the conditions described above.

    Raises:
        `ValueError`: when ``proteins`` contain domains without probabilities.

    Example:
        In the following example, we check the significance of three domains
        (*A*, *B* and *C*) on BGC membership for a training set containing
        7 proteins:

        >>> data = [
        ...     Protein("prot1", _, [
        ...         Domain("A", _, _, _, _, _, probability=1),
        ...         Domain("B", _, _, _, _, _, probability=1),
        ...     ]),
        ...     Protein("prot2", _, [
        ...         Domain("A", _, _, _, _, _, probability=1),
        ...         Domain("B", _, _, _, _, _, probability=1),
        ...     ]),
        ...     Protein("prot3", _, [
        ...         Domain("A", _, _, _, _, _, probability=1),
        ...         Domain("B", _, _, _, _, _, probability=1),
        ...     ]),
        ...     Protein("prot4", _, [Domain("A", _, _, _, _, _, probability=1)]),
        ...     Protein("prot5", _, [Domain("A", _, _, _, _, _, probability=1)]),
        ...     Protein("prot6", _, [
        ...         Domain("C", _, _, _, _, _, probability=0),
        ...         Domain("B", _, _, _, _, _, probability=0),
        ...     ]),
        ...     Protein("prot7", _, [Domain("C", _, _, _, _, _, probability=0)]),
        ... ]
        >>> sorted(fisher_significance(data).items())
        [('A', 0.071...), ('B', 0.999...), ('C', 0.071...)]

        Since *A* and *C* only appear in BGC and non BGC proteins respectively,
        the p-value for a two-tailed Fisher Exact Test is under 5%, while *B*,
        which appears in half of the BGC proteins and in half of the non-BGC
        proteins, is not significant with regards to the fisher test.

        It's also possible to get the uncorrected values by giving `None`
        instead of a correction method:

        >>> sorted(fisher_significance(data, correction_method=None).items())
        [('A', 0.047...), ('B', 0.999...), ('C', 0.047...)]

    """
    # set of all proteins, +proteins grouped by features
    proteins_ = set(), set()  # type: ignore
    features_ = collections.defaultdict(set), collections.defaultdict(set)  # type: ignore

    # collect proteins / features for all data tables
    for protein in proteins:
        for domain in protein.domains:
            if domain.probability is None:
                raise ValueError("domain is missing a BGC probability")
            in_bgc = domain.probability > 0.5
            proteins_[in_bgc].add(protein.id)
            features_[in_bgc][domain.name].add(protein.id)

    # make the contigency table for each feature
    significance = {}
    for feature in set(features_[False]).union(features_[True]):
        pvalue = fisher.pvalue(
            len(features_[True][feature]),  # with feature, in BGC
            len(proteins_[True]) - len(features_[True][feature]),  # without feature, in BGC
            len(features_[False][feature]),  # with feature, not in BGC
            len(proteins_[False]) - len(features_[False][feature]),  # without feature, not in BGC
        )
        significance[feature] = pvalue.two_tail

    # perform multiple test correction if needed
    if correction_method is not None:
        significance = significance_correction(significance, correction_method)

    return significance
