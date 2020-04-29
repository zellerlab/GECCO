"""Features extraction and sanitization utilities for `~gecco.crf.ClusterCRF`.
"""

import collections
import math
import numbers
import operator
import typing
import warnings
from collections.abc import Iterable
from itertools import zip_longest
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy
import pandas


def truncate(
    data: "pandas.DataFrame",
    length: int,
    label_column: str = "BGC",
    group_column: str = "protein_id",
) -> "pandas.DataFrame":
    """Truncate ``data`` while preserving positively labeled regions.

    Arguments:
        data (`~pandas.DataFrame`): The data table to truncate.
        length (`int`): The target length for the new data table.
        label_column (`str`): The name of the column containing labels in the
            data table.
        group_column (`str`): The name of the column to use to group the data.
            Rows of the same group are kept together or excluded together.

    """
    data0 = data[data[label_column] == 0]
    data1 = data[data[label_column] == 1]
    data0 = [df for _, df in data0.groupby(group_column, sort=False)]
    trunc_len = (len(data0) - 2 * length) // 2

    # try:
    data0_trunc = pandas.concat(data0[trunc_len:-trunc_len])
    # except ValueError as err:
    #    data_trunc = pandas.concat(data0)

    return data0_trunc.append(data1).sort_index().reset_index()


def extract_group_features(
    table: "pandas.DataFrame",
    feature_columns: List[str],
    weight_columns: List[str],
    group_column: List[str],
    label_column: Optional[str] = None,
) -> Tuple[List[Dict[str, float]], Optional[List[str]]]:
    """Extract features from ``table`` on a group level.

    Extraction is done respecting the ``feature_columns`` and ``weight_columns``
    arguments on each group obtained with ``group_column``.

    This function is mostly used to group on a *protein* level, but it can
    potentially be used as well to group on a larger subunit, for instance
    on properly labeled contiguous ORFs to mimick a BGC.

    Arguments:
        table (~pandas.DataFrame): The dataframe to process.
        feature_columns (iterable of `str`): The name of the columns containing
            the features in the table.
        weight_columns (iterable of `str`): The name of the columns containing
            the weights in the table.
        label_column (`str`, optional): The name of the column containing the class
            labels, or `None`.
        overlap (`int`): The size of the sliding window; features for the
            row at index *n* will be extracted from rows located at index in
            *[n-overlap; n+overlap]*.

    Returns:
        `tuple`: a couple of `list`, where the first list contains a
        feature dictionary for each group, and the second the list of
        class labels, or `None` if ``label_column`` was `None`.

    Example:
        >>> data = pandas.DataFrame(
        ...     columns=["protein_id", "domain", "weight"],
        ...     data=[
        ...         ["prot1", "domainA", 0.5],
        ...         ["prot1", "domainB", 1.0],
        ...         ["prot2", "domainC", 0.8],
        ...     ],
        ... )
        >>> extract_group_features(data, ["domain"], ["weight"], "protein_id")
        ([{'domainA': 0.5, 'domainB': 1.0}, {'domainC': 0.8}], None)

    """
    # NOTE(@althonos):
    # There is a lot to unpack here, but basically we are doing some convoluted
    # stuff in order to avoid calling `table.groupby` directly, which takes
    # ages otherwise. To do so, we work at the `numpy.ndarray` level, extracting
    # the columns and doing the grouping at the same time we extract features.
    # The code also does some branching to avoid extracting labels if it is
    # not necessary.

    # first, let's extract all the groups, and create a table that maps each
    # group to an index: this lets us store the features of each groups in an
    # array instead of a dictionary, sparing us the list conversion at the end
    groups = list(table[group_column].unique())
    group_indexer = {group: i for i, group in enumerate(groups)}
    X: List[Dict[str, float]] = [dict() for _ in groups]

    # then let's extract all the arrays we need from the input dataframe
    # so that we can discard it and only work with numpy arrays
    group_array = table[group_column].array
    feature_weight_arrays = [
        (table[feat].array, table[weight].array)
        for feat, weight in zip(feature_columns, weight_columns)
    ]

    if label_column is None:
        # iterate on all the rows of the initial dataframe, and for
        # each row extract all `feature:weight` couples and put them in
        # the right group, which is determined by the group column.
        for row_index, group_key in enumerate(group_array):
            x = X[group_indexer[group_key]]
            for feat_array, weight_array in feature_weight_arrays:
                row_feat = feat_array[row_index]
                row_weight = weight_array[row_index]
                x[row_feat] = max(x.get(row_feat, 0), row_weight)

        return X, None

    else:
        # if we are extracting class labels, we also need the label array
        label_array = table[label_column].array.astype(str)
        Y: List[Optional[str]] = [None for _ in groups]

        # then, same as before, but with some extract code at the end
        # to extract the class label, and check a group does not contain
        # mismatching labels
        for row_index, group_key in enumerate(group_array):
            group_index = group_indexer[group_key]
            x, y = X[group_index], Y[group_index]
            for feat_array, weight_array in feature_weight_arrays:
                row_feat = feat_array[row_index]
                row_weight = weight_array[row_index]
                x[row_feat] = max(x.get(row_feat, 0), row_weight)

            label_row = label_array[row_index]
            if y is not None and y != label_row:
                warnings.warn(
                    "Feature group contains mixed class label. A random label will be selected."
                )
            else:
                Y[group_index] = label_row

        return X, typing.cast(List[str], Y)


def extract_overlapping_features(
    table: "pandas.DataFrame",
    feature_columns: List[str],
    weight_columns: List[str],
    overlap: int = 1,
    label_column: Optional[str] = None,
) -> Tuple[List[Dict[str, float]], Optional[List[str]]]:
    """Extract features from ``table`` using a sliding window.

    The extraction is done using the ``feature_columns`` and ``weight_columns``
    columns to extract features on each overlapping window.

    Arguments:
        table (~pandas.DataFrame): The dataframe to process.
        feature_columns (iterable of `str`): The name of the columns containing
            the features in the table.
        weight_columns (iterable of `str`): The name of the columns containing
            the weights in the table.
        label_column (`str`, optional): The name of the column containing the
            class labels, or `None`.
        overlap (`int`): The size of the sliding window; features for the
            row at index *n* will be extracted from rows located at index in
            *[n-overlap; n+overlap]*.

    Returns:
        `tuple`: a couple of `list`, where the first list contains a
        feature dictionary for each row, and the second the list of
        class labels, or `None` if ``label_column`` was `None`.

    Example:
        >>> data = pandas.DataFrame(
        ...     columns=["protein_id", "domain", "weight"],
        ...     data=[
        ...         ["prot1", "A", 0.5],
        ...         ["prot1", "B", 1.0],
        ...         ["prot2", "C", 0.8],
        ...     ],
        ... )
        >>> extract_overlapping_features(data, ["domain"], ["weight"])[0]
        [{'A': 0.5, 'B': 1.0}, {'A': 0.5, 'B': 1.0, 'C': 0.8}, {'B': 1.0, 'C': 0.8}]

    """
    # create a feature list for each slice of the data
    X: List[Dict[str, float]] = [dict() for _ in range(len(table))]
    for idx in range(len(table)):
        # get the indices of the sliding window
        start_idx = max(idx - overlap, 0)
        end_idx = min(idx + overlap + 1, len(table))
        # process the features
        for feat_col, weight_col in zip(feature_columns, weight_columns):
            features = table[feat_col].array[start_idx:end_idx]
            weights = table[weight_col].array[start_idx:end_idx]
            for feat, weight in zip(features, weights):
                X[idx][feat] = max(X[idx].get(feat, 0), weight)
    # Only return Y if requested
    if label_column is None:
        return X, None
    return X, list(table[label_column].array.astype(str))


def extract_single_features(
    table: "pandas.DataFrame",
    feature_columns: List[str],
    weight_columns: List[str],
    label_column: Optional[str] = None,
) -> Tuple[List[Dict[str, float]], Optional[List[str]]]:
    """Extract features from ``table`` on a row level.

    The extraction is done using the ``feature_columns`` and ``weight_columns``
    columns to extract features on each row.

    Arguments:
        table (~pandas.DataFrame): The dataframe to process.
        feature_columns (iterable of `str`): The name of the columns containing
            the features in the table.
        weight_columns (iterable of `str`): The name of the columns containing
            the weights in the table.
        label_column (`str`, optional): The name of the column containing the
            class labels, or `None`.

    Returns:
        `tuple`: a couple of `list`, where the first list contains a
        feature dictionary for each row, and the second the list of
        class labels, or `None` if ``label_column`` was `None`.

    Example:
        >>> data = pandas.DataFrame(
        ...     columns=["protein_id", "domain", "weight"],
        ...     data=[
        ...         ["prot1", "domainA", 0.5],
        ...         ["prot1", "domainB", 1.0],
        ...         ["prot2", "domainC", 0.8],
        ...     ],
        ... )
        >>> extract_single_features(data, ["domain"], ["weight"])
        ([{'domainA': 0.5}, {'domainB': 1.0}, {'domainC': 0.8}], None)

    """
    # extract weights without iterating on all rows by zipping together
    # the appropriate columns and inserting them in the right location
    X: List[Dict[str, float]] = [dict() for _ in range(len(table))]
    for feat_col, weight_col in zip(feature_columns, weight_columns):
        features, weights = table[feat_col].array, table[weight_col].array
        for index, (feature, weight) in enumerate(zip(features, weights)):
            X[index][feature] = weight
    # return Y only if a label column is given
    if label_column is None:
        return X, None
    return X, list(table[label_column].array.astype(str))
