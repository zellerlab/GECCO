import math
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.model_selection import PredefinedSplit
from sklearn_crfsuite import CRF
from itertools import zip_longest
from orion.cross_validation import LotoSplit, n_folds, n_folds_partial, StratifiedSplit
from orion.preprocessing import extract_features, flatten, truncate
from orion.preprocessing import extract_overlapping_features, extract_protein_features


# CLASS
class ClusterCRF(object):

    def __init__(self, Y_col=None,
        feature_cols=[], weight_cols=[], group_col="protein_id", feature_type="single",
        algorithm="lbsgf", overlap=2, **kwargs):

        self.Y_col = Y_col
        self.features = feature_cols
        self.weights = weight_cols
        self.groups = group_col
        self.feature_type = feature_type
        self.overlap = overlap

        self.alg = algorithm
        self.model = CRF(
            algorithm = self.alg,
            all_possible_transitions = True,
            all_possible_states = True,
            **kwargs)

    def fit(self, X=None, Y=None, data=None):
        if (X is not None) and (Y is not None):
            self.model.fit(X, Y)
        elif data is not None:
            samples = [self._extract_features(s) for s in data]
            X = np.array([x for x, _ in samples])
            Y = np.array([y for _, y in samples])
            self.model.fit(X, Y)

    def predict_marginals(self, data=None, X=None):
        if X is not None:
            return self.model.predict_marginals(X)
        elif data is not None:
            samples = [self._extract_features(s) for s in data]
            X = np.array([x for x, _ in samples])
            marginal_probs = self.model.predict_marginals(X)
            marginal_probs = np.concatenate(
                np.array([np.array(_) for _ in marginal_probs]))
            cluster_probs = np.array([d["1"] for d in [s for s in marginal_probs]])

            if self.feature_type == "group":
                data = [self._distinct(s) for s in data]

            result_df = pd.concat(data)
            result_df = result_df.assign(p_pred=cluster_probs)

            return result_df

    def cv(self, data, k=10, threads=1, e_filter=1, trunc=None, strat_col=None):

        if strat_col:
            types = [s[strat_col].values[0].split(",") for s in data]
            cv_split = StratifiedSplit(types, n_splits=k)
        else:
            folds = n_folds(len(data), n=k)
            cv_split = PredefinedSplit(folds)

        results = Parallel(n_jobs=threads)(
            delayed(self._single_fold_cv)(data, train_idx, test_idx, e_filter=e_filter,
                trunc=trunc) for train_idx, test_idx in cv_split.split())

        return results

    def loto_cv(self, data, type_col, threads=1, e_filter=1, trunc=None):

        types = [s[type_col].values[0].split(",") for s in data]
        cv_split = LotoSplit(types)

        results = Parallel(n_jobs=threads)(
            delayed(self._single_fold_cv)(data, train_idx, test_idx,
                round_id=typ, e_filter=e_filter, trunc=trunc)
                for train_idx, test_idx, typ in cv_split.split())

        return results

    def partial_cv(self, data, n_train, n_val, k=10, threads=1, e_filter=1,
            trunc=None):

        folds = n_folds_partial(n_train, n_val, n=k)
        cv_split = PredefinedSplit(folds)

        results = Parallel(n_jobs=threads)(
            delayed(self._single_fold_cv)(data, train_idx, test_idx, e_filter=e_filter,
                trunc=trunc) for train_idx, test_idx in cv_split.split())

        return results

    def _single_fold_cv(self, data, train_idx, test_idx, round_id=None, e_filter=1,
            trunc=None):
        """Performs a single CV round with the given train_idx and test_idx
        """

        train_data = [data[i].reset_index() for i in train_idx]

        if truncate:
            train_data = [truncate(df, trunc, Y_col=self.Y_col, grouping=self.groups)
                for df in train_data]

        train_samples = [self._extract_features(s) for s in train_data]
        X_train = np.array([x for x, _ in train_samples])
        Y_train = np.array([y for _, y in train_samples])

        test_data = [data[i].reset_index() for i in test_idx]
        test_data = [t[t["i_Evalue"] < e_filter].reset_index(drop=True)
            for t in test_data]
        test_samples = [self._extract_features(s) for s in test_data]

        X_test = np.array([x for x, _ in test_samples])
        Y_test = np.array([y for _, y in test_samples])

        self.fit(X=X_train, Y=Y_train)

        marginal_probs = self.model.predict_marginals(X_test)
        marginal_probs = np.concatenate(np.array([np.array(_) for _ in marginal_probs]))
        cluster_probs = np.array([d["1"] for d in [s for s in marginal_probs]])

        Y_pred = np.array([1 if p > 0.5 else 0 for p in cluster_probs])
        Y_real = np.array([int(i) for i in np.concatenate(Y_test)])

        if self.feature_type == "group":
            test_data = [self._distinct(s) for s in test_data]

        result_df = pd.concat(test_data)
        result_df = result_df.assign(
            p_pred = cluster_probs,
            y_pred = Y_pred,
            cv_round = round_id)

        return result_df

    def _extract_features(self, sample):

        if self.feature_type == "single":
            return extract_features(sample, self.Y_col,
                feature_col=self.features,
                weight_col=self.weights)

        if self.feature_type == "overlap":
            return extract_overlapping_features(sample, self.Y_col,
                feature_col=self.features,
                weight_col=self.weights,
                overlap=self.overlap)

        if self.feature_type == "group":
            return extract_protein_features(sample, self.Y_col,
                feature_col=self.features,
                weight_col=self.weights,
                prot_col=self.groups)

    def _distinct(self, df):
        return  df.groupby(self.groups, sort=False).first().reset_index()
