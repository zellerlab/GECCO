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


# CLASS
class ClusterCRF(object):
    """A wrapper for `sklearn_crfsuite.CRF` taking `~pandas.DataFrame` inputs.

    `ClusterCRF` enables prediction and cross-validation for dataframes. This
    is handy to use with feature tables obtained from `~gecco.hmmer.HMMER`.
    """

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

        if feature_type == "single":
            self._extract_features = self._extract_single_features
        elif feature_type == "group":
            self._extract_features = self._extract_group_features
        elif feature_type == "overlap":
            self._extract_features = self._extract_overlapping_features
        else:
            raise ValueError(f"invalid feature type: {feature_type!r}")

    def fit(self, X=None, Y=None, data=None):
        """Fits the model to the given data.

        If X and Y are defined, it fits the model on the corresponding vectors.
        If data is defined, it takes a dataframe to extract features and labels for
        fitting.
        """
        if X is not None and Y is not None and data is not None:
            raise ValueError("either X and Y or data must be given, not both")

        if data is not None:
            samples = [self._extract_features(s) for s in data]
            X = numpy.array([x for x, _ in samples])
            Y = numpy.array([y for _, y in samples])
            self.model.fit(X, Y)
        else:
            self.model.fit(X, Y)

    def predict_marginals(self, data=None, X=None):
        """Predicts marginals for your data.

        If X (a feature vector) is defined, it outputs a probability vector.
        If data (a dataframe) is defined, it outputs the same dataframe with the
        probability vector concatenated to it as the column p_pred.
        """
        if X is not None:
            return self.model.predict_marginals(X)
        elif data is not None:
            # Extract features to dict (CRFSuite format)
            samples = [self._extract_features(s, X_only=True) for s in data]
            X = numpy.array([x for x, _ in samples])
            marginal_probs = self.model.predict_marginals(X)
            # Extract cluster (1) probs
            cluster_probs = []
            for sample in marginal_probs:
                try:
                    sample_probs = numpy.array([d["1"] for d in [_ for _ in sample]])
                except KeyError:
                    warnings.warn(
                        "Cluster probabilities of test set were found to be zero. This indicates that there might be something wrong with your input data.", Warning
                    )
                    sample_probs = numpy.array([0 for d in [_ for _ in sample]])
                cluster_probs.append(sample_probs)
            # Merge probs vector with the input dataframe. This is tricky if we are
            # dealing with protein features as length of vector does not fit to dataframe
            # To deal with this, we merge by protein_id
            # --> HOWEVER: this requires the protein IDs to be unique among all samples
            if self.feature_type == "group":
                result_df = [self._merge(df, p_pred=p) for df, p in zip(data, cluster_probs)]
                result_df = pandas.concat(result_df)
            else:
                result_df = pandas.concat(data)
                result_df = result_df.assign(p_pred=numpy.concatenate(cluster_probs))

            return result_df

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
        """Performs a single CV round with the given train_idx and test_idx.
        """

        # Extract the CV fold from the complete data using the provided indices
        test_data = [data[i].reset_index() for i in test_idx]
        train_data = [data[i].reset_index() for i in train_idx]
        if trunc is not None:
            # Truncate training set from both sides to desired length
            train_data = [
                truncate(df, trunc, Y_col=self.Y_col, grouping=self.groups)
                for df in train_data
            ]

        # Extract features to CRFSuite format
        with multiprocessing.pool.Pool(threads) as pool:
            # only keep useful columns to reduce the total amount of data
            # passed to the other processes
            columns = self.features + self.weights + [self.groups, self.Y_col]
            col_filter = operator.itemgetter(columns)

            train_samples = pool.map(self._extract_features, map(col_filter, train_data))
            X_train = numpy.array([x for x, _ in train_samples])
            Y_train = numpy.array([y for _, y in train_samples])

            test_samples = pool.map(self._extract_features, map(col_filter, test_data))
            X_test = numpy.array([x for x, _ in test_samples])
            Y_test = numpy.array([y for _, y in test_samples])

        # Fit the model
        self.fit(X=X_train, Y=Y_train)

        # Extract cluster (1) probabilities from marginals
        marginal_probs = self.model.predict_marginals(X_test)
        cluster_probs = [
            numpy.array([d.get("1", 0) for d in sample])
            for sample in marginal_probs
        ]

        # check if any sample has P(1) == 0
        if any(not prob.all() for prob in cluster_probs):
            warnings.warn("Cluster probabilities of test set were found to be zero. Something may be wrong with your input data.")

        # Merge probs vector with the input dataframe. This is tricky if we are dealing
        # with protein features as length of vector does not fit to dataframe
        # To deal with this, we merge by protein_id
        # --> HOWEVER: this requires the protein IDs to be unique within each sample
        if self.feature_type == "group":
            result_df = pandas.concat([self._merge(df, p_pred=p) for df, p in zip(test_data, cluster_probs)])
        else:
            result_df = pandas.concat(test_data).assign(p_pred=numpy.concatenate(cluster_probs))
        return result_df.assign(cv_round=round_id)

    def _merge(self, df, **cols):
        unidf = pandas.DataFrame(cols)
        unidf[self.groups] = df[self.groups].unique()
        return df.merge(unidf)

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
                X[-1].update(zip(df[feat_col], df[weight_col]))

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
