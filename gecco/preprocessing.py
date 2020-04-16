import math
import numbers
import typing
from collections import Iterable
from itertools import zip_longest
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.base import TransformerMixin, BaseEstimator


class SafeOneHotEncoder(BaseEstimator, TransformerMixin):
    """
    Onehot encodes samples with differing categorical composition
    given the full set of categories.
    """
    def __init__(self, all_categories):
        self.le = LabelEncoder()
        self.ohe = OneHotEncoder()
        self.ohe.fit(self.le.fit_transform(all_categories).reshape(-1,1))

    def transform(self, X):
        return self.ohe.transform(self.le.transform(X).reshape(-1,1)).toarray()


def flatten(l):
    """Flattens list of arbitraty depth to generator."""
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def truncate(df, length, label_column="BGC", group_column="protein_id"):

    df0 = df[df[label_column] == 0]
    df1 = df[df[label_column] == 1]
    df0 = [df for _, df in df0.groupby(group_column, sort=False)]
    trunc_len = (len(df0) - 2 * length) // 2

    try:
        df_trunc = pandas.concat(df0[trunc_len : -trunc_len])
    except ValueError as err:
        df_trunc = pandas.concat(df0)

    return df_trunc.append(df1).sort_index().reset_index()


def safe_encode(X, feature_set=None, encoding="onehot"):
    """Sefely onehot encodes samples with different categorical composition."""
    if isinstance(feature_set, Iterable):
        categories = list(feature_set)
    else:
        categories = list(set(flatten(X)))
    if encoding == "onehot":
        enc = SafeOneHotEncoder(categories)
    elif encoding == "label":
        enc = LabelEncoder().fit(categories)
    X_enc = np.array([enc.transform(_.ravel()) for _ in X])
    if isinstance(feature_set, Iterable):
        return X_enc
    else:
        return X_enc, categories


def extract_group_features(
        table: pandas.DataFrame,
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
        >>> extract_group_features(data, ["domain"], ["weight"], ["protein_id"])
        ([{"domainA": 0.5, "domainB": 1.0}, {"domainC": 0.8}], None)

    """
    # create a feature list for each group (i.e. protein) without
    # iterating on each row
    X, Y = [], []
    for prot_id, df in table.groupby(group_column, sort=False):
        X.append({})
        for feat_col, weight_col in zip(feature_columns, weight_columns):
            features, weights = df[feat_col].values, df[weight_col].values
            for feat, weight in zip(features, weights):
                X[-1][feat] = max(X[-1].get(feat, 0), weight)
        if label_column is not None:
            Y.append(str(df[label_column].values[0]))
    # only return Y if the class column was given
    return X, None if label_column is None else Y



def extract_overlapping_features(
    table: pandas.DataFrame,
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
        ...         ["prot1", "domainA", 0.5],
        ...         ["prot1", "domainB", 1.0],
        ...         ["prot2", "domainC", 0.8],
        ...     ],
        ... )
        >>> extract_overlapping_features(data, ["domain"], ["weight"])
        (
            [
                {"domainA": 0.5, "domainB": 1.0},
                {"domainA": 0.5, "domainB": 1.0, "domainC": 0.8},
                {"domainB": 1.0, "domainC": 0.8}
            ],
            None
        )

    """
    # create a feature list for each slice of the data
    X = [dict() for _ in range(len(table))]
    for idx in range(len(table)):
        # get the indices of the sliding window
        start_idx = max(idx-overlap, 0)
        end_idx = min(idx+overlap+1, len(table))
        # process the features
        for feat_col, weight_col in zip(feature_columns, weight_columns):
            features = table[feat_col].values[start_idx:end_idx]
            weights = table[weight_col].values[start_idx:end_idx]
            for feat, weight in zip(features, weights):
                X[idx][feat] = max(X[idx].get(feat, 0), weight)
    # Only return Y if requested
    return X, None if label_column is None else table[label_column].values.astype(str)


def extract_single_features(
    table: pandas.DataFrame,
    feature_columns: List[str],
    weight_columns: List[str],
    label_column: Optional[str] = None
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
        ([{"domainA": 0.5}, {"domainB": 1.0}, {"domainC": 0.8}], None)

    """
    # extract weights without iterating on all rows by zipping together
    # the appropriate columns and inserting them in the right location
    X = [dict() for _ in range(len(table))]
    for feat_col, weight_col in zip(feature_columns, weight_columns):
        features, weights = table[feat_col].values, table[weight_col].values
        for index, (feature, weight) in enumerate(zip(features, weights)):
            X[index][feature] = weight
    # return Y only if a label column is given
    return X, None if label_column is None else table[label_column].values.astype(str)
