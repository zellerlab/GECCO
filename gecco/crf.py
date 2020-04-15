import csv
import functools
import operator
import math
import multiprocessing.pool
import numbers
import random
import warnings
from typing import List, Optional

import numpy
import pandas
import tqdm
from sklearn.model_selection import PredefinedSplit
from sklearn_crfsuite import CRF

from .cross_validation import LotoSplit, n_folds, n_folds_partial, StratifiedSplit
from .preprocessing import flatten, truncate


class ClusterCRF(object):
    """A wrapper for `sklearn_crfsuite.CRF` taking `~pandas.DataFrame` inputs.

    `ClusterCRF` enables prediction and cross-validation for dataframes. This
    is handy to use with feature tables obtained from `~gecco.hmmer.HMMER`. It
    supports arbitrary column names that can be changed on initialisation.
    """

    _EXTRACT_FEATURES_METHOD = {
        "single": "_extract_single_features",
        "group": "_extract_group_features",
        "overlap": "_extract_overlapping_features"
    }

    def __init__(
            self,
            Y_col: Optional[str] = None,
            feature_cols: Optional[List[str]] = None,
            weight_cols=None,
            group_col="protein_id",
            feature_type="single",
            algorithm="lbsgf",
            overlap=2,
            **kwargs
    ) -> None:
        """Create a new `ClusterCRF` instance.

        Arguments:
            Y_col (`str`): The name of the column containing class labels. Must
                be given if the model is going to be trained, but not needed
                when only making predictions.
            feature_cols (list of `str`): The name of the column(s) with
                categorical features.
            weight_cols (list of `str`): The name of the column(s) with
                weights for categorical features. *These are applied locally
                and don't correspond to the actual weights the model learns.*
                See also the `~sklearn_crfsuite.CRFSuite` documentation.
            group_col (str): In case of `feature_type = "group"`,  defines the
                grouping column to use.
            feature_type (str): Defines how features should be extracted. The
                following values are accepted:

                - ``single``: features are extracted on a domain/row level
                - ``overlap``: features are extracted in overlapping windows
                - ``group``: features are extracted in groupings determined
                  by a column in the data frame. *This is most useful when
                  dealing with proteins, but can handle arbitrary grouping
                  levels*.

            algorithm (str): The optimization algorithm for the model. See
                https://sklearn-crfsuite.readthedocs.io/en/latest/api.html
                for available values.
            overlap (int): In case of `feature_type = "overlap"`, defines the
                window size to use.

        Any additional keyword argument is passed as-is to the internal
        `~sklearn_crfsuite.CRF` constructor.
        """
        if feature_type not in self._EXTRACT_FEATURES_METHOD:
            raise ValueError(f"unexpected feature type: {feature_type!r}")

        self.Y_col = Y_col
        self.features = feature_cols or []
        self.weights = weight_cols or []
        self.groups = group_col
        self.feature_type = feature_type
        self.overlap = overlap
        self.algorithm = algorithm
        self.model = CRF(
            algorithm = algorithm,
            all_possible_transitions = True,
            all_possible_states = True,
            **kwargs
        )

    def fit(self, data, threads=None):
        """Fits the model to the given data.

        Arguments:
            data (iterable of `~pandas.DataFrame`): An iterable of data frames
                to use to fit the model. If must contain the following columns
                (as defined on `ClusterCRF` initialisation): *weight_cols*,
                *feature_cols* and *Y_col*, as well as *feature_cols* if the
                feature extraction is ``group``.
            threads (`int`, optional): The number of threads to use when
                extracting the features from the input data.

        """
        X, Y = self._extract_features(data, threads=threads)
        self.model.fit(X, Y)

    def predict_marginals(self, data, threads=None):
        """Predicts marginals for your data.

        Arguments:
            data (iterable of `~pandas.DataFrame`): An iterable of data frames
                to use to fit the model. If must contain the following columns
                (as defined on `ClusterCRF` initialisation): *weight_cols* and
                *feature_cols*, as well as *feature_cols* if the
                feature extraction is ``group``.
            threads (`int`, optional): The number of threads to use when
                extracting the features from the input data.

        """
        # convert data to `CRFSuite` format
        X, _ = self._extract_features(data, threads=threads, X_only=True)

        # Extract cluster (1) probabilities from marginal
        marginal_probs = self.model.predict_marginals(X)
        cluster_probs = [
            numpy.array([d.get("1", 0) for d in sample])
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

    def cv(self, data, strat_col=None, k=10, threads=1, trunc=None):
        """Runs k-fold cross-validation using a stratification column.

        Arguments:
            data (`~pandas.DataFrame`): A domain annotation table.
            k (int): The number of cross-validation fold to perform.
            threads (int): The number of threads to use.
            trunc (int, optional): The maximum number of rows to use in the
                training data, or to `None` to use everything.
            strat_col (str, optional): The name of the column to use to split
                the data, or `None` to perform a predefined split.

        Returns:
            `list` of `~pandas.DataFrame`: The list containing one table of
            results for each cross-validation fold.

        Todo:
            * Fix multiprocessing but within `sklearn` to launch each fold in
              a separate `~multiprocessing.pool.ThreadPool`.
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
                threads=threads
            )
            for i, (train_idx, test_idx) in enumerate(pbar)
        ]

    def loto_cv(self, data, strat_col, threads=1, trunc=None):
        """Run LOTO cross-validation using a stratification column.

        Arguments:
            data (`~pandas.DataFrame`): A domain annotation table.
            strat_col (str): The name of the column to use to split the data.
            threads (int): The number of threads to use.
            trunc (int, optional): The maximum number of rows to use in the
                training data, or to `None` to use everything.

        Returns:
            `list` of `~pandas.DataFrame`: The list containing one table for
            each cross-validation fold.

        Todo:
            * Fix multiprocessing but within `sklearn` to launch each fold in
              a separate `~multiprocessing.pool.ThreadPool`.
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

    def _single_fold_cv(self, data, train_idx, test_idx, round_id=None, trunc=None, threads=None):
        """Performs a single CV round with the given indices.
        """
        # Extract the fold from the complete data using the provided indices
        train_data = [data[i].reset_index() for i in train_idx]
        if trunc is not None:
            # Truncate training set from both sides to desired length
            train_data = [
                truncate(df, trunc, Y_col=self.Y_col, grouping=self.groups)
                for df in train_data
            ]

        # Fit the model
        self.fit(train_data, threads=threads)

        # Predict marginals on test data and return predictions
        test_data = [data[i].reset_index() for i in test_idx]
        marginals = self.predict_marginals(data=test_data, threads=threads)
        return marginals.assign(cv_round=round_id)

    # --- Feature extraction -------------------------------------------------

    def _extract_features(self, data, threads=None, X_only=False):
        """Convert a data list to `CRF`-compatible wrappers.

        Arguments:
            data (`list` of `~pandas.DataFrame`): A list of samples to
                extract features and class labels from.
            threads (`int`, optional): The number of parallel threads to
                launch to extract features.
            X_only (`bool`, optional): If `True`, only return the features,
                and ignore class labels.

        Warning:
            This method spawns a `multiprocessing.pool.Pool` in the background
            to extract features in parallel.
        """
        # Filter the columns to reduce the amount of data passed to the
        # different processes
        columns = self.features + [self.groups]
        columns.extend(filter(lambda w: isinstance(w, str), self.weights))
        if not X_only:
            columns.append(self.Y_col)
        col_filter = operator.itemgetter(columns)
        # Extract features to CRFSuite format
        _extract_features = functools.partial(
            getattr(self, self._EXTRACT_FEATURES_METHOD[self.feature_type]),
            X_only=X_only,
        )
        with multiprocessing.pool.Pool(threads) as pool:
            samples = pool.map(_extract_features, map(col_filter, data))
            X = numpy.array([x for x, _ in samples])
            Y = numpy.array([y for _, y in samples])
        # Only return Y if requested
        return X, None if X_only else Y

    def _extract_single_features(self, table, X_only=False):
        """
        Prepares class labels Y and features from a table
        given
        table:  a table with pfam domains and class lables
        Y_col:  the column encoding the class labels (e.g. 1 (BGC) and 0 (non-BGC))
        if Y_col == None, only X is returned --> for prediction cases
        """
        X = []
        for _, row in table.iterrows():
            feat_dict = dict()
            feat_dict = self._make_feature_dict(row, feat_dict)
            X.append(feat_dict)
        if not X_only:
            Y = numpy.array(table[self.Y_col].astype(str))
            return X, Y
        else:
            return X, None

    def _extract_group_features(self, table, X_only=False):
        """Extract features from ``table`` on a group level.

        The extraction is done respecting the ``group_col``, ``feature_cols``
        and ``weight_cols`` arguments given on the `ClusterCRF` initialisation.

        This function is mostly used to group on a *protein* level, but it can
        potentially be used as well to group on a larger subunit, for instance
        on properly labeled contiguous ORFs to mimick a BGC.

        Arguments:
            table (~pandas.DataFrame): The dataframe to process.
            X_only (bool): If `True`, prevents extraction of class labels
                from the groups.

        Returns:
            `tuple`: a couple of `list`, where the first list contains a
            feature dictionary for each group, and the second the list of
            class labels, or `None` if ``class_column`` was none.

        Example:
            >>> data = pandas.DataFrame(
            ...     columns=["protein_id", "domain", "weight"],
            ...     data=[
            ...         ["prot1", "domainA", 0.5],
            ...         ["prot1", "domainB", 1.0],
            ...         ["prot2", "domainC", 0.8],
            ...     ],
            ... )
            >>> crf = ClusterCRF(
            ...     feature_cols=["domain"],
            ...     group_col="protein_id",
            ...     feature_type="group",
            ... )
            >>> crf._extract_group_features(data, X_only=True)
            ([{"domainA": 0.5, "domainB": 1.0}, {"domainC": 0.8}], None)

        """

        # if we have only one feature to extract, we can sort the data outside
        # of the loop to greatly decrease the processing within the loop:
        # we sort and drop duplicates so that each feature appears only once
        # per group, with its maximum weight
        if len(self.features) == 1:
            table = (
                table
                    .sort_values(by=[self.groups] + self.features + self.weights)
                    .drop_duplicates([self.groups] + self.features, keep="last")
            )

        # create a feature list for each group (i.e. protein)
        X, Y = [], []
        for prot_id, df in table.groupby(self.groups, sort=False):
            X.append({})
            if not X_only:
                Y.append(str(df[self.Y_col].values[0]))
            for feat_col, weight_col in zip(self.features, self.weights):
                # if we couldn't preprocess because we have more than one
                # column here, we need to sort the values so that the
                # biggest weight comes last
                if len(self.features) > 1:
                    df = df.sort_values(by=[feat_col, weight_col])
                # we can now create the feature dictionary, using a Python
                # implementation detail: when adding several time the same
                # key to a dictionary with `update`, then only the last one
                # is kept: because we sorted the last key is the biggest
                # weight
                X[-1].update(zip(df[feat_col], df.get(weight_col, weight_col)))

        # only return Y if the class column was given
        return X, None if X_only else Y

    def _extract_overlapping_features(self, table, X_only=False):
        """
        Prepares class labels Y and features from a table
        given
        table:  a table with pfam domains and class lables
        Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
        if Y_col == None, only X is returned --> for prediction cases
        """
        X = []
        for idx, _ in table.iterrows():
            wind = table.iloc[idx - self.overlap : idx + self.overlap + 1]
            feat_dict = dict()
            for _, row in wind.iterrows():
                feat_dict = self._make_feature_dict(row, feat_dict)
            X.append(feat_dict)
        if not X_only:
            Y = numpy.array(table[self.Y_col].astype(str))
            return X, Y
        else:
            return X, None

    def _make_feature_dict(self, row, feat_dict=None):
        """
        Constructs a dict with key:value pairs from
        row: input row or dict
        self.features: either name of feature or name of column in row/dict
        self.weights: either numerical weight or name of column in row/dict
        """
        feat_dict = feat_dict or {}
        for f, w in zip(self.features, self.weights):
            if isinstance(w, numbers.Number):
                key, val = row[f], w
            else:
                key, val = row.get(f, f), row[w]
            feat_dict[key] = max(val, feat_dict.get(key, val))
        return feat_dict

    # --- Utils --------------------------------------------------------------

    def _merge(self, df, **cols):
        unidf = pandas.DataFrame(cols)
        unidf[self.groups] = df[self.groups].unique()
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
