import csv
import functools
import operator
import math
import multiprocessing.pool
import numbers
import random
import warnings
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union

import numpy
import pandas
import tqdm
from sklearn.model_selection import PredefinedSplit
from sklearn_crfsuite import CRF

from . import preprocessing
from .cross_validation import LotoSplit, n_folds, n_folds_partial, StratifiedSplit


class ClusterCRF(object):
    """A wrapper for `sklearn_crfsuite.CRF` taking `~pandas.DataFrame` inputs.

    `ClusterCRF` enables prediction and cross-validation for dataframes. This
    is handy to use with feature tables obtained from `~gecco.hmmer.HMMER`. It
    supports arbitrary column names that can be changed on initialisation.
    """

    def __init__(
        self,
        feature_columns: Union[str, List[str]],
        weight_columns: Union[str, List[str]],
        group_column: str = "protein_id",
        label_column: str = "BGC",
        feature_type: str = "single",
        algorithm: str = "lbsgf",
        overlap: int = 2,
        **kwargs,
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
            group_column (`str`): If ``feature_type`` is *group*`, defines
                the column to use to group rows together. *This can be, for
                instance, the column containing protein IDs to extract features
                on the same gene together*.
            label_column (`str`): The name of the column containing class labels
                (the **Y** column). Not actually used when making predictions.
            feature_type (str): Defines how features should be extracted. The
                following values are accepted:

                - ``single``: features are extracted on a domain/row level
                - ``overlap``: features are extracted in overlapping windows
                - ``group``: features are extracted in groupings determined
                  by a column in the data frame. *This is most useful when
                  dealing with proteins, but can handle arbitrary grouping
                  levels*.

            algorithm (`str`): The optimization algorithm for the model. See
                https://sklearn-crfsuite.readthedocs.io/en/latest/api.html
                for available values.
            overlap (`int`): In case of `feature_type = "overlap"`, defines
                the sliding window size to use. The resulting window width is
                ``2*overlap+1``.

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
        self.model = CRF(
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
        jobs: Optional[int]=None,
    ) -> None:
        """Fits the model to the given data.

        Arguments:
            data (iterable of `~pandas.DataFrame`): An iterable of data frames
                to use to fit the model. If must contain the following columns
                (as defined on `ClusterCRF` initialisation): ``weight_columns``,
                ``feature_columns`` and ``label_column``, as well as
                ``group_column`` if the feature extraction is ``group``.
            jobs (`int`, optional): The number of jobs to use to extract the
                features from the input data.

        """
        X, Y = self._extract_features(data, jobs=jobs)
        self.model.fit(X, Y)

    def predict_marginals(
        self,
        data: Iterable[pandas.DataFrame],
        jobs: Optional[int]=None
    ) -> pandas.DataFrame:
        """Predicts marginals for your data.

        Arguments:
            data (iterable of `~pandas.DataFrame`): An iterable of data frames
                to use to fit the model. If must contain the following columns
                (as defined on `ClusterCRF` initialisation): *weight_cols* and
                *feature_cols*, as well as *feature_cols* if the
                feature extraction is ``group``.
            jobs (`int`, optional): The number of jobs to use to extract the
                features from the input data.

        """
        # convert data to `CRFSuite` format
        X, _ = self._extract_features(data, jobs=jobs, X_only=True)

        # Extract cluster (1) probabilities from marginal
        marginal_probs = self.model.predict_marginals(X)
        cluster_probs = [
            numpy.array([data.get("1", 0) for data in sample])
            for sample in marginal_probs
        ]

        # check if any sample has P(1) == 0
        if any(not prob.all() for prob in cluster_probs):
            warnings.warn(
                """
                Cluster probabilities of test set were found to be zero.
                Something may be wrong with your input data.
                """
            )

        # Merge probs vector with the input dataframe. This is tricky if we are
        # dealing with protein features as length of vector does not fit to dataframe
        # To deal with this, we merge by protein_id
        # --> HOWEVER: this requires the protein IDs to be unique among all samples
        if self.feature_type == "group":
            results = [self._merge(df, p_pred=p) for df, p in zip(data, cluster_probs)]
            return pandas.concat(results)
        else:
            return pandas.concat(data).assign(p_pred=numpy.concatenate(cluster_probs))

    # --- Cross-validation ---------------------------------------------------

    def cv(
        self,
        data: List[pandas.DataFrame],
        strat_col: Optional[str] = None,
        k: int = 10,
        jobs: Optional[int] = None,
        trunc: Optional[int]=None
    ) -> List[pandas.DataFrame]:
        """Runs k-fold cross-validation using a stratification column.

        Arguments:
            data (`list` of `~pandas.DataFrame`): A list of domain annotation
                table, corresponding to different samples.
            k (`int`): The number of cross-validation folds to perform.
            jobs (`int`): The number of jobs to use to extract the features from
                the input data within each fold.
            trunc (`int`, optional): The maximum number of rows to use in the
                training data, or to `None` to use everything.
            strat_col (`str`, optional): The name of the column to use to split
                the data, or `None` to perform a predefined split.

        Returns:
            `list` of `~pandas.DataFrame`: The list containing one table of
            results for each cross-validation fold.

        Todo:
            * Make progress bar configurable.
        """
        if strat_col is not None:
            types = [s[strat_col].values[0].split(",") for s in data]
            cv_split = StratifiedSplit(types, n_splits=k)
        else:
            folds = n_folds(len(data), n=k)
            cv_split = PredefinedSplit(folds)

        # Not running in parallel because sklearn has issues managing the
        # temporary files in multithreaded mode
        pbar = tqdm.tqdm(cv_split.split(), total=k, leave=False)
        return [
            self._single_fold_cv(
                data,
                train_idx,
                test_idx,
                round_id=f"fold{i}",
                trunc=trunc,
                jobs=jobs,
            )
            for i, (train_idx, test_idx) in enumerate(pbar)
        ]

    def loto_cv(
        self,
        data: List[pandas.DataFrame],
        strat_col: str,
        jobs: Optional[int]=None,
        trunc: Optional[int]=None
    ) -> List[pandas.DataFrame]:
        """Run LOTO cross-validation using a stratification column.

        Arguments:
            data (`list` of `~pandas.DataFrame`): A list of domain annotation
                table, corresponding to different samples.
            strat_col (`str`): The name of the column to use to split the data.
            jobs (`int`): The number of jobs to use to extract the features from
                the input data within each fold.
            trunc (`int`, optional): The maximum number of rows to use in the
                training data, or to `None` to use everything.

        Returns:
            `list` of `~pandas.DataFrame`: The list containing one table for
            each cross-validation fold.

        Todo:
            * Make progress bar configurable.
        """
        labels = [s[strat_col].values[0].split(",") for s in data]
        cv_split = LotoSplit(labels)

        # Not running in parallel because sklearn has issues managing the
        # temporary files in multithreaded mode
        pbar = tqdm.tqdm(list(cv_split.split()), leave=False)
        return [
            self._single_fold_cv(data, train_idx, test_idx, label, trunc)
            for train_idx, test_idx, label in pbar
        ]

    def _single_fold_cv(
        self,
        data: List[pandas.DataFrame],
        train_idx: numpy.ndarray,
        test_idx: numpy.ndarray,
        round_id: Optional[str] = None,
        trunc: Optional[int]=None,
        jobs: Optional[int]=None
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
        self.fit(train_data, jobs=jobs)

        # Predict marginals on test data and return predictions
        test_data = [data[i].reset_index() for i in test_idx]
        marginals = self.predict_marginals(data=test_data, jobs=jobs)
        if round_id is not None:
            marginals["cv_round"] = round_id
        return marginals

    # --- Feature extraction -------------------------------------------------

    def _currify_extract_function(
        self,
        X_only: bool = False
    ) -> Callable[
        [pandas.DataFrame],
        Tuple[List[Dict[str, float]], Optional[List[str]]]
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
                preprocessing.extract_group_features,  # type: ignore
                group_column=self.group_column,
            )
        elif self.feature_type == "single":
            extract = preprocessing.extract_single_features  # type: ignore
        elif self.feature_type == "overlap":
            extract = functools.partial(
                preprocessing.extract_overlapping_features,  # type: ignore
                overlap=self.overlap
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
        data: List[pandas.DataFrame],
        jobs: Optional[int] = None,
        X_only: bool = False
    ) -> Tuple[List[Dict[str, float]], Optional[List[str]]]:
        """Convert a data list to `CRF`-compatible wrappers.

        Arguments:
            data (`list` of `~pandas.DataFrame`): A list of samples to
                extract features and class labels from.
            jobs (`int`, optional): The number of jobs to use to extract the
                features from the input data.
            X_only (`bool`, optional): If `True`, only return the features,
                and ignore class labels.

        Warning:
            This method spawns a `multiprocessing.pool.Pool` in the background
            to extract features in parallel.
        """
        # Filter the columns to reduce the amount of data passed to the
        # different processes
        columns = self.feature_columns + self.weight_columns + [self.group_column]
        if not X_only:
            columns.append(self.label_column)
        col_filter = operator.itemgetter(columns)
        # Extract features to CRFSuite format
        _extract = self._currify_extract_function(X_only=X_only)
        with multiprocessing.pool.Pool(jobs) as pool:
            samples = pool.map(_extract, map(col_filter, data))
            X = numpy.array([x for x, _ in samples])
            Y = numpy.array([y for _, y in samples])
        # Only return Y if requested
        return X, None if X_only else Y

    # --- Utils --------------------------------------------------------------

    def _merge(self, df, **cols):
        unidf = pandas.DataFrame(cols)
        unidf[self.group_column] = df[self.group_column].unique()
        return df.merge(unidf)

    def save_weights(self, basename: str) -> None:
        with open(f"{basename}.trans.tsv", "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["from", "to", "weight"])
            for labels, weight in self.model.transition_features_.items():
                writer.writerow([*labels, weight])
        with open(f"{basename}.state.tsv", "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["attr", "label", "weight"])
            for attrs, weight in self.model.state_features_.items():
                writer.writerow([*attrs, weight])
