"""BGC prediction using a conditional random field.
"""

import csv
import functools
import hashlib
import itertools
import io
import math
import numbers
import operator
import os
import pickle
import random
import textwrap
import typing
import warnings
from multiprocessing.pool import Pool
from typing import Any, Callable, Dict, FrozenSet, Iterable, List, Optional, Tuple, Type, Union

import numpy
import pandas
import tqdm
import pkg_resources
import sklearn_crfsuite
import sklearn.model_selection
import sklearn.preprocessing

from . import preprocessing
from .cv import LeaveOneGroupOut
from .select import fisher_significance


class ClusterCRF(object):
    """A wrapper for `sklearn_crfsuite.CRF` taking `~pandas.DataFrame` inputs.

    `ClusterCRF` enables prediction and cross-validation for dataframes. This
    is handy to use with feature tables obtained from `~gecco.hmmer.HMMER`. It
    supports arbitrary column names that can be changed on initialisation.

    Extraction of features can be done in different modes, which change the
    level at which features are grouped. *Note that each model is specific for
    a ``feature_type`` value, since changing the extraction method will lead
    to different training data for the model*. The following methods are
    supported:

    ``single``
        Features are extracted on a row level, which corresponds to single
        domains most of the time.
    ``overlap``
        Features are extracted in overlapping windows, allowing to group
        together contiguous domains.
    ``group``
        Features are extracted in groupings determined
        by a column in the data frame. *This is most useful when
        dealing with proteins, but can handle arbitrary grouping
        levels*.

    Warning:
        This class attempts to parallelize some operations, namely extracting
        features of several samples. By default, it uses processes, which can
        achieve true parallelism contrary to threads, but it's not always
        possible to launch new processes, for instance when the `ClusterCRF` is
        invoked within a daemon process. If this happen, try giving a different
        class to the ``pool_factory`` argument of the constructor, like
        `~multiprocessing.pool.ThreadPool`.

    """

    @classmethod
    def trained(cls, model_path: Optional[str] = None) -> "ClusterCRF":
        """Create a new pre-trained `ClusterCRF` instance from a model file.

        Arguments:
            model_path (`str`, optional): The path to the model directory
                obtained with the ``gecco train`` command. If `None` given,
                use the embedded model.

        Returns:
            `~gecco.crf.ClusterCRF`: A CRF model that can be used to perform
            predictions without training first.

        Raises:
            `ValueError`: If the model data does not match its hash.

        """
        if model_path is not None:
            pkl_path = os.path.join(model_path, "model.pkl")
            md5_path = os.path.join(model_path, "model.pkl.md5")
        else:
            pkl_path = pkg_resources.resource_filename(__name__, "model.pkl")
            md5_path = pkg_resources.resource_filename(__name__, "model.pkl.md5")
        with open(md5_path) as sig:
            signature = sig.read().strip()

        hasher = hashlib.md5()
        with open(pkl_path, "rb") as bin:
            read = functools.partial(bin.read, io.DEFAULT_BUFFER_SIZE)
            for chunk in iter(read, b""):
                hasher.update(typing.cast(bytes, chunk))
        if hasher.hexdigest().upper() != signature.upper():
            raise ValueError("MD5 hash of model data does not match signature")

        with open(pkl_path, "rb") as bin:
            return pickle.load(bin)  # type: ignore

    def __init__(
        self,
        feature_columns: Union[str, List[str]],
        weight_columns: Union[str, List[str]],
        group_column: str = "protein_id",
        label_column: str = "BGC",
        feature_type: str = "single",
        algorithm: str = "lbfgs",
        overlap: int = 2,
        pool_factory: Union[Type[Pool], Callable[[Optional[int]], Pool]] = Pool,
        **kwargs: Dict[str, object],
    ) -> None:
        """Create a new `ClusterCRF` instance.

        Arguments:
            feature_columns (`str`, or iterable of `str`): The name(s) of the
                column(s) with categorical features.
            weight_columns (`str`, or iterable of `str`): The name(s) of the
                column(s) with weights for categorical features. *These are
                applied locally and don't correspond to the actual weights
                the model learns.* See also the `~sklearn_crfsuite.CRFSuite`
                documentation.
            group_column (`str`): If ``feature_type`` is *group*, defines
                the column to use to group rows together. *This can be, for
                instance, the column containing protein IDs to extract features
                on the same gene together*.
            label_column (`str`): The name of the column containing class labels
                (the **Y** column). Not actually used when making predictions.
            feature_type (`str`): Defines how features should be extracted.
                Should be one *single*, *group*, or *overlap*.
            algorithm (`str`): The optimization algorithm for the model. See
                https://sklearn-crfsuite.readthedocs.io/en/latest/api.html
                for available values.
            overlap (`int`): In case of `feature_type = "overlap"`, defines
                the sliding window size to use. The resulting window width is
                ``2*overlap+1``.
            pool_factory (`multiprocessing.pool.Pool` subclass, or callable):
                The factory to use to create a new pool instance for methods
                that can perform operations in parallel. *It is called with
                a single argument which is the number of workers to create,
                or* `None` *to create a much workers as there are CPUs.*

        Any additional keyword argument is passed as-is to the internal
        `~sklearn_crfsuite.CRF` constructor.

        Raises:
            `ValueError`: if ``feature_type`` has an invalid value.
            `TypeError`: if one of the ``*_columns`` argument is not iterable.

        """
        if feature_type not in {"single", "overlap", "group"}:
            raise ValueError(f"unexpected feature type: {feature_type!r}")

        self.feature_type: str = feature_type
        self.label_column: str = label_column
        self.group_column: str = group_column
        self.overlap: int = overlap
        self.algorithm = algorithm
        self.pool_factory = pool_factory
        self.significant_features: Dict[str, FrozenSet[str]]  = {}
        self.model = sklearn_crfsuite.CRF(
            algorithm=algorithm,
            all_possible_transitions=True,
            all_possible_states=True,
            **kwargs,
        )

        if isinstance(feature_columns, str):
            self.feature_columns = [feature_columns]
        else:
            self.feature_columns = list(feature_columns)

        if isinstance(weight_columns, str):
            self.weight_columns = [weight_columns]
        else:
            self.weight_columns = list(weight_columns)

    def fit(
        self,
        data: Iterable[pandas.DataFrame],
        trunc: Optional[int] = None,
        select: Optional[float] = None,
        *,
        jobs: Optional[int] = None,
    ) -> None:
        """Fits the model to the given data.

        Arguments:
            data (iterable of `~pandas.DataFrame`): An iterable of data frames
                to use to fit the model. If must contain the following columns
                (as defined on `ClusterCRF` initialisation): ``weight_columns``,
                ``feature_columns`` and ``label_column``, as well as
                ``group_column`` if the feature extraction method is *group*.
            trunc (`int`, optional): A number of maximum rows to truncate the
                data groups to.
            select (`float`, optional): The fraction of most significant
                features to keep, or `None` to disable feature selection. For
                instance, a value of :math:`0.1` means the 10% most
                significant features will be used for training and prediction.

        Keyword Arguments:
            jobs (`int`, optional): The number of jobs to use to extract the
                features from the input data. If `None` given, use as many jobs
                as there are CPUs.

        """
        # Collect data to avoid iteration issues
        data = list(data)

        # Truncate data if requested
        if trunc is not None:
            data = [
                preprocessing.truncate(df, trunc, self.label_column, self.group_column)
                for df in data
            ]

        # Find most significant domains in the training data
        # QUESTION(@althonos): should it be before or after data truncation ?
        if select is not None:
            if select <= 0 or select > 1:
                raise ValueError(f"invalid value for select: {select}")
            for feat in self.feature_columns:
                sig = fisher_significance(data, feat, self.label_column, self.group_column)
                sorted_sig = sorted(sig, key=sig.get)[:int(select*len(sig))]
                self.significant_features[feat] = frozenset(sorted_sig)
            for column, features in self.significant_features.items():
                data = [ df[ df[column].isin(features) ] for df in data ]

        # Do the actual training
        X, Y = self._extract_features(data, jobs=jobs)
        self.model.fit(X, Y)

    def predict_marginals(
        self, data: Iterable[pandas.DataFrame], *, jobs: Optional[int] = None,
    ) -> pandas.DataFrame:
        """Predicts marginals for the input data.

        Arguments:
            data (iterable of `~pandas.DataFrame`): An iterable of data frames
                to predict classes for. If must contain the following columns
                (as defined on `ClusterCRF` initialisation): ``weight_columns``
                and ``feature_columns``, as well as ``group_column`` if the
                feature extraction method is *group*.

        Keyword Arguments:
            jobs (`int`, optional): The number of jobs to use to extract the
                features from the input data. If `None` given, use as many jobs
                as there are CPUs.

        """
        # Convert data to `CRFSuite` format
        X, _ = self._extract_features(data, X_only=True, jobs=jobs)

        # Remove non-significant features from the extracted bunch
        if self.significant_features:
            sf = set(itertools.chain(*self.significant_features.values()))
            X = [[{k: row[k] for k in row.keys() & sf} for row in x] for x in X]

        # Extract cluster (1) probabilities from predicted marginals
        marginal_probs = self.model.predict_marginals(X)
        cluster_probs = [
            numpy.array([data.get("1", 0) for data in sample])
            for sample in marginal_probs
        ]

        # Check if any sample has P(1) == 0
        if any(not prob.all() for prob in cluster_probs):
            warnings.warn(
                textwrap.dedent(
                    """
                    Cluster probabilities of test set were found to be zero.
                    Something may be wrong with your input data.
                    """
                )
            )

        # Merge probs vector with the input dataframe. This is tricky if we are
        # dealing with protein features as length of vector does not fit to dataframe
        # To deal with this, we merge by protein_id
        # --> HOWEVER: this requires the protein IDs to be unique among all samples
        if self.feature_type == "group":
            results = [self._merge(df, p_pred=p) for df, p in zip(data, cluster_probs)]  # type: ignore
            return pandas.concat(results)
        else:
            return pandas.concat(data).assign(p_pred=numpy.concatenate(cluster_probs))

    # --- Cross-validation ---------------------------------------------------

    def cv(
        self,
        data: List[pandas.DataFrame],
        strat_col: Optional[str] = None,
        k: int = 10,
        trunc: Optional[int] = None,
        select: Optional[float] = None,
        *,
        jobs: Optional[int] = None,
    ) -> List[pandas.DataFrame]:
        """Runs k-fold cross-validation, possibly with a stratification column.

        Arguments:
            data (`list` of `~pandas.DataFrame`): A list of domain annotation
                table, corresponding to different samples.
            strat_col (`str`, optional): The name of the column to use to
                stratify the cross-validation folds or `None` to simply split
                the input data in ``k`` folds.
            k (`int`): The number of cross-validation folds to perform.
            trunc (`int`, optional): The maximum number of rows to use in the
                training data, or to `None` to use everything.
            select (`float`, optional): The fraction of most significant
                features to keep, or `None` to disable feature selection. For
                instance, a value of :math:`0.1` means the 10% most
                significant features will be used for training and prediction.

        Keyword Arguments:
            jobs (`int`): The number of jobs to use to extract the features from
                the input data within each fold. If `None` given, use as many
                jobs as there are CPUs.

        Returns:
            `list` of `~pandas.DataFrame`: A list containing one table of
            results (predicted on the testing data) for each cross-validation
            fold.

        Todo:
            * Make progress bar configurable.
        """
        if strat_col is not None:
            cross_validator = sklearn.model_selection.StratifiedKFold(k)
            # we need to extract the BGC types (given in `s[strat_col]`) and
            # to linearize them (using "Mixed" as a type if a sequence has
            # more than one BGC type)
            strat_labels = [s[strat_col].values[0].split(";") for s in data]
            y = ["Mixed" if len(label) > 1 else label[0] for label in strat_labels]
            splits = list(cross_validator.split(data, y=y))
        else:
            cross_validator = sklearn.model_selection.KFold(k)
            splits = list(cross_validator.split(data))

        return [
            self._single_fold_cv(
                data,
                train_idx=trn,
                test_idx=tst,
                round_id=f"fold{i}",
                trunc=trunc,
                select=select,
                jobs=jobs,
            )
            for i, (trn, tst) in enumerate(tqdm.tqdm(splits, leave=False))
        ]

    def loto_cv(
        self,
        data: List[pandas.DataFrame],
        strat_col: str,
        trunc: Optional[int] = None,
        select: Optional[float] = None,
        *,
        jobs: Optional[int] = None,
    ) -> List[pandas.DataFrame]:
        """Run LOTO cross-validation using a stratification column.

        Arguments:
            data (`list` of `~pandas.DataFrame`): A list of domain annotation
                table, corresponding to different samples.
            strat_col (`str`): The name of the column to use to stratify
                the cross-validation folds.
            trunc (`int`, optional): The maximum number of rows to use in the
                training data, or to `None` to use everything.
            select (`float`, optional): The fraction of most significant
                features to keep, or `None` to disable feature selection. For
                instance, a value of :math:`0.1` means the 10% most
                significant features will be used for training and prediction.

        Keyword Arguments:
            jobs (`int`): The number of jobs to use to extract the features from
                the input data within each fold. If `None` given, use as many
                jobs as there are CPUs.

        Returns:
            `list` of `~pandas.DataFrame`: A list containing one table of
            results (predicted on the testing data) for each cross-validation
            fold.

        Todo:
            * Make progress bar configurable.
        """
        cross_validator = LeaveOneGroupOut()
        strat_labels = [tuple(s[strat_col].values[0].split(";")) for s in data]
        splits = list(cross_validator.split(data, groups=strat_labels))

        return [
            self._single_fold_cv(
                data,
                train_idx=trn,
                test_idx=tst,
                round_id=f"fold{i}",
                trunc=trunc,
                select=select,
                jobs=jobs,
            )
            for i, (trn, tst) in enumerate(tqdm.tqdm(splits, leave=False))
        ]

    def _single_fold_cv(
        self,
        data: List[pandas.DataFrame],
        train_idx: numpy.ndarray,
        test_idx: numpy.ndarray,
        round_id: Optional[str] = None,
        trunc: Optional[int] = None,
        select: Optional[float] = None,
        *,
        jobs: Optional[int] = None,
    ) -> pandas.DataFrame:
        """Performs a single CV round with the given indices.
        """
        # Extract the fold from the complete data using the provided indices
        train_data = [data[i].reset_index() for i in train_idx]
        if trunc is not None:
            # Truncate training set from both sides to desired length
            train_data = [
                preprocessing.truncate(df, trunc, self.label_column, self.group_column)
                for df in train_data
            ]

        # Fit the model
        self.fit(train_data, jobs=jobs, select=select)

        # Predict marginals on test data and return predictions
        test_data = [data[i].reset_index() for i in test_idx]
        marginals = self.predict_marginals(data=test_data, jobs=jobs)
        if round_id is not None:
            marginals["cv_round"] = round_id
        return marginals

    # --- Feature extraction -------------------------------------------------

    def _currify_extract_function(
        self, X_only: bool = False
    ) -> Callable[
        [pandas.DataFrame], Tuple[List[Dict[str, float]], Optional[List[str]]]
    ]:
        """Currify a feature extraction function from `gecco.preprocessing`.

        Dependending on the `feature_type` given on the `CRFSuite`
        initialisation, different extraction functions can be used. Since
        they do not have the same signature, they can't be swapped with each
        other easily. This method currifies these functions so that they
        all take a list of `~pandas.DataFrame` as single positional argument.
        """
        if self.feature_type == "group":
            extract = functools.partial(
                preprocessing.extract_group_features, group_column=self.group_column,
            )
        elif self.feature_type == "single":
            extract = functools.partial(preprocessing.extract_single_features)
        elif self.feature_type == "overlap":
            extract = functools.partial(
                preprocessing.extract_overlapping_features, overlap=self.overlap,
            )
        else:
            raise ValueError(f"invalid feature type: {self.feature_type!r}")
        return functools.partial(
            extract,
            feature_columns=self.feature_columns,
            weight_columns=self.weight_columns,
            label_column=None if X_only else self.label_column,
        )

    def _extract_features(
        self,
        data: Iterable[pandas.DataFrame],
        X_only: bool = False,
        *,
        jobs: Optional[int] = None,
    ) -> Tuple[List[List[Dict[str, float]]], Optional[List[List[str]]]]:
        """Convert a data list to `CRF`-compatible wrappers.

        Arguments:
            data (`list` of `~pandas.DataFrame`): A list of samples to
                extract features and class labels from.
            X_only (`bool`, optional): If `True`, only return the features,
                and ignore class labels.

        Keyword Arguments:
            jobs (`int`, optional): The number of jobs to use to extract the
                features from the input data.

        """
        # Filter the columns to reduce the amount of data passed to the
        # different processes
        columns: List[str] = self.feature_columns + self.weight_columns
        if self.feature_type == "group":
            columns.append(self.group_column)
        if not X_only:
            columns.append(self.label_column)
        col_filter = operator.itemgetter(columns)
        # Extract features to CRFSuite format
        _extract = self._currify_extract_function(X_only=X_only)
        with self.pool_factory(jobs) as pool:
            samples = pool.map(_extract, map(col_filter, data))
            X, Y = [x for x, _ in samples], [y for _, y in samples]
        # Only return Y if requested
        return X, None if X_only else typing.cast(List[List[str]], Y)

    # --- Utils --------------------------------------------------------------

    def _merge(self, df: pandas.DataFrame, **cols: List[Any]) -> pandas.DataFrame:
        unidf = pandas.DataFrame(cols)
        unidf[self.group_column] = df[self.group_column].unique()
        return df.merge(unidf)

    def save_weights(self, directory: str) -> None:
        with open(os.path.join(directory, "model.trans.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["from", "to", "weight"])
            for labels, weight in self.model.transition_features_.items():
                writer.writerow([*labels, weight])
        with open(os.path.join(directory, "model.trans.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["attr", "label", "weight"])
            for attrs, weight in self.model.state_features_.items():
                writer.writerow([*attrs, weight])
