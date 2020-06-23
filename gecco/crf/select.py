"""Feature selection using Fisher's Exact Test.
"""

import collections
import typing
from typing import Dict, Iterable, Mapping, Optional

import fisher
import pandas
import numpy
from statsmodels.stats.multitest import fdrcorrection

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
    data: Iterable["~pandas.DataFrame"],
    feature_column: str = "domain",
    label_column: str = "BGC",
    protein_column: str = "protein_id",
    correction_method: Optional[str] = "indep",
) -> Dict[str, float]:
    r"""Estimate the significance of each feature on the label prediction.

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
        data (`~pandas.DataFrame`): The table containing the entirety of the
            training set.
        feature_column (`str`): The name of the column containing feature
            names in the table. *Only a single feature column is supported,
            but this function can be called multiple times on different
            feature columns if needed*.
        label_column (`str`): The column which contains the class label, which
            should be either *0* or *1* depending on whether the protein is in
            a BGC.
        protein_column (`str`): The name of the column containing protein
            identifiers, since contigency counting occurs at the protein level
            independently of the feature extraction method used for the CRF.
        correction_method (`str`, optional): The name of the multiple test
            correction method to use when computing significance, or `None` to
            skip correction. See `statsmodels.stats.multitest.fdrcorrection`
            for allowed values.

    Returns:
        dict: A dictionary which to each feature associates the p-value of
        the two-tail Fisher Exact Test in the conditions described above.

    Example:
        In the following example, we check the significance of three domains
        (*A*, *B* and *C*) on BGC membership for a training set containing
        7 proteins:

        >>> data = pandas.DataFrame(
        ...     columns=["protein_id", "domain", "BGC"],
        ...     data=[
        ...         ["prot1", "A", 1], ["prot1", "B", 1],
        ...         ["prot2", "A", 1], ["prot2", "B", 1],
        ...         ["prot3", "A", 1], ["prot3", "B", 1],
        ...         ["prot4", "A", 1],
        ...         ["prot5", "A", 1],
        ...         ["prot6", "C", 0], ["prot6", "B", 0],
        ...         ["prot7", "C", 0],
        ...     ],
        ... )
        >>> sorted(fisher_significance([data]).items())
        [('A', 0.071...), ('B', 0.999...), ('C', 0.071...)]

        Since *A* and *C* only appear in BGC and non BGC proteins respectively,
        the p-value for a two-tailed Fisher Exact Test is under 5%, while *B*,
        which appears in half of the BGC proteins and in half of the non-BGC
        proteins, is not significant with regards to the fisher test.

        It's also possible to get the uncorrected values by giving `None`
        instead of a correction method:

        >>> sorted(fisher_significance([data], correction_method=None).items())
        [('A', 0.047...), ('B', 0.999...), ('C', 0.047...)]

    """
    # set of all proteins, +proteins grouped by features
    proteins = set(), set()
    features = collections.defaultdict(set), collections.defaultdict(set)

    # collect proteins / features for all data tables
    for df in data:
        # directly work on numpy arrays to avoid using `groupby`
        array_prot = df[protein_column].values
        array_feat = df[feature_column].values
        array_bgc = df[label_column].values
        # for each feature, extract whether they are in a BGC or not, and store
        # the identifier of the protein they are part of
        for prot, feature, in_bgc in zip(array_prot, array_feat, array_bgc):
            proteins[in_bgc].add(prot)
            features[in_bgc][feature].add(prot)

    # make the contigency table for each feature
    significance = {}
    for feature in set(features[0]).union(features[1]):
        pvalue = fisher.pvalue(
            len(features[1][feature]),  # proteins with feature, in BGC
            len(proteins[1]) - len(features[1][feature]), # proteins without feature, in BGC
            len(features[0][feature]), # proteins with feature, not in BGC
            len(proteins[0]) - len(features[0][feature]),
        )
        significance[feature] = pvalue.two_tail

    # perform multiple test correction if needed
    if correction_method is not None:
        significance = significance_correction(significance, correction_method)

    return significance