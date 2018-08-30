import math
import warnings
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.model_selection import PredefinedSplit
from sklearn_crfsuite import CRF
from itertools import zip_longest
from gecco.cross_validation import LotoSplit, n_folds, n_folds_partial, StratifiedSplit
from gecco.preprocessing import extract_features, flatten, truncate
from gecco.preprocessing import extract_overlapping_features, extract_protein_features


# CLASS
class ClusterCRF(object):
    """
    ClusterCRF is a wrapper around CRFSuite and enables predition and cross validation
    for dataframes. This is handy for sequence prediction in e.g. annotated genomic
    sequences.
    Initiated with the columns defining the respective properties in the dataframe:
        Y_col: column with class labels
        feature_cols: column(s) with categorical features
        weight_cols: column(s) with 'weights' for categorical features. These are applied
            locally and don't correspont to the actual weights the model learns. You can read up on how this is used inside CRFSuite on the CRFSuite website.
        group_col: in case of feature_type = 'group', this defines the grouping column
        feature_type: defines how features should be extracted:
            single: Features are extractd on a domain (row) level
            overlap: Features are extracted in overlapping windows
            group: Features are extracted in groupings determined by a column in the data
            frame. This is most useful when dealing with proteins, but can handle
            arbitrary grouping levels
        algorithm: the optimization algorithm for the model
            (again, check the CRFSuite website)
        overlap: in case of feature_type='overlap', by how much the windows should
            overlap
        **kwargs: other parameters you want to pass to the CRF model.
    """

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
        """
        Fits the model.
        If X and Y are defined, it fits the model on the corresponding vectors.
        If data is defined, it takes toe dataframe to extract features and labels for
        fitting
        """
        if (X is not None) and (Y is not None):
            self.model.fit(X, Y)
        elif data is not None:
            samples = [self._extract_features(s) for s in data]
            X = np.array([x for x, _ in samples])
            Y = np.array([y for _, y in samples])
            self.model.fit(X, Y)

    def predict_marginals(self, data=None, X=None):
        """
        Predicts marginals for your data.
        If X (a feature vector) is defined, it outputs a probability vector.
        If data (a dataframe) is defined, it outputs the same dataframe with the
        probability vector concatenated to it as the column p_pred.
        """
        if X is not None:
            return self.model.predict_marginals(X)
        elif data is not None:
            # Extract features to dict (CRFSuite format)
            samples = [self._extract_features(s, X_only=True) for s in data]
            X = np.array([x for x, _ in samples])
            marginal_probs = self.model.predict_marginals(X)
            # Extract cluster (1) probs
            marginal_probs = np.concatenate(
                np.array([np.array(_) for _ in marginal_probs]))
            cluster_probs = np.array([d["1"] for d in [s for s in marginal_probs]])

            # Merge probs vector with the input dataframe. This is tricky if we are
            # dealing with protein features as length of vector does not fit to dataframe
            # To deal with this, we merge by protein_id
            # --> HOWEVER: this requires the protein IDs to be unique among all samples
            if self.feature_type == "group":
                groups = np.concatenate([df[self.groups].unique() for df in data])
                result_df = pd.concat(data)
                result_df = self._merge(result_df, groups, p_pred=cluster_probs)
            else:
                result_df = pd.concat(data)
                result_df = result_df.assign(p_pred=cluster_probs)

            return result_df

    def cv(self, data, k=10, threads=1, e_filter=1, trunc=None, strat_col=None):
        """Runs stratified k-fold CV using k splits and a stratification column"""

        if strat_col:
            types = [s[strat_col].values[0].split(",") for s in data]
            cv_split = StratifiedSplit(types, n_splits=k)
        else:
            folds = n_folds(len(data), n=k)
            cv_split = PredefinedSplit(folds)

        # Run one job per split and collects the result in a list
        results = Parallel(n_jobs=threads)(
            delayed(self._single_fold_cv)(data, train_idx, test_idx, e_filter=e_filter,
                trunc=trunc) for train_idx, test_idx in cv_split.split())
        return results

    def loto_cv(self, data, type_col, threads=1, e_filter=1, trunc=None):
        """Runs LOTO CV based on a column defining the type of the sample (type_col)"""

        types = [s[type_col].values[0].split(",") for s in data]
        cv_split = LotoSplit(types)

        # Run one job per split and collects the result in a list
        results = Parallel(n_jobs=threads)(
            delayed(self._single_fold_cv)(data, train_idx, test_idx,
                round_id=typ, e_filter=e_filter, trunc=trunc)
                for train_idx, test_idx, typ in cv_split.split())
        return results

    def partial_cv(self, data, n_train, n_val, k=10, threads=1, e_filter=1,
            trunc=None):
        """
        Runs partial CV. A training set is only used for training and the
        validation set is used for k-fold CV
        """

        folds = n_folds_partial(n_train, n_val, n=k)
        cv_split = PredefinedSplit(folds)

        # Run one job per split and collects the result in a list
        results = Parallel(n_jobs=threads)(
            delayed(self._single_fold_cv)(data, train_idx, test_idx, e_filter=e_filter,
                trunc=trunc) for train_idx, test_idx in cv_split.split())
        return results

    def _single_fold_cv(self, data, train_idx, test_idx, round_id=None, e_filter=1,
            trunc=None):
        """Performs a single CV round with the given train_idx and test_idx
        """

        train_data = [data[i].reset_index() for i in train_idx]

        if trunc:
            # Truncate training set from both sides to desired length
            train_data = [truncate(df, trunc, Y_col=self.Y_col, grouping=self.groups)
                for df in train_data]

        # Extract features to CRFSuite format
        train_samples = [self._extract_features(s) for s in train_data]
        X_train = np.array([x for x, _ in train_samples])
        Y_train = np.array([y for _, y in train_samples])

        test_data = [data[i].reset_index() for i in test_idx]
        # Optional filter for E-value
        test_data = [t[t["i_Evalue"] < e_filter].reset_index(drop=True)
            for t in test_data]
        test_samples = [self._extract_features(s) for s in test_data]

        X_test = np.array([x for x, _ in test_samples])
        Y_test = np.array([y for _, y in test_samples])

        self.fit(X=X_train, Y=Y_train)

        # Extract cluster (1) probabilities from marginals
        marginal_probs = self.model.predict_marginals(X_test)
        marginal_probs = np.concatenate(np.array([np.array(_) for _ in marginal_probs]))
        try:
            cluster_probs = np.array([d["1"] for d in [s for s in marginal_probs]])
        except KeyError:
            warnings.warn(
                "Cluster probabilities of test set were found to be zero. This indicates that there might be something wrong with your input data.", Warning
            )
            cluster_probs = np.array([0 for d in [s for s in marginal_probs]])

        # Merge probs vector with the input dataframe. This is tricky if we are dealing
        # with protein features as length of vector does not fit to dataframe
        # To deal with this, we merge by protein_id
        # --> HOWEVER: this requires the protein IDs to be unique among all samples
        if self.feature_type == "group":
            groups = np.concatenate([df[self.groups].unique() for df in test_data])
            result_df = pd.concat(test_data).assign(cv_round=round_id)
            result_df = self._merge(result_df, groups, p_pred=cluster_probs)
        else:
            result_df = (pd.concat(test_data)
                .assign(
                    p_pred = cluster_probs,
                    cv_round = round_id
                )
            )
        return result_df

    def _extract_features(self, sample, X_only=False):
        """
        Chooses extraction function based on feature type.
        single: Features are extractd on a domain (row) level
        overlap: Features are extracted in overlapping windows
        group: Features are extracted in groupings determined by a column in the data
        frame. This is most useful when dealing with proteins, bu can handle arbitrary
        grouping levels
        """
        # Little hacky this, but... meh...
        if X_only:
            Y_col = None
        else:
            Y_col = self.Y_col

        if self.feature_type == "single":
            return self._extract_features(sample, Y_col)

        if self.feature_type == "overlap":
            return self._extract_overlapping_features(sample, Y_col)

        if self.feature_type == "group":
            return self._extract_protein_features(sample, Y_col)

    def _merge(self, df, groups, **cols):
        unidf = pd.DataFrame(cols)
        unidf[self.groups] = groups
        return df.merge(unidf)

    def _extract_features(self, table, Y_col=None):
        """
        Prepares class labels Y and features from a table
        given
        table:  a table with pfam domains and class lables
        Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
        """
        X = []
        for _, row in table.iterrows():
            feat_dict = dict()
            feat_dict = self._make_feature_dict(row, feat_dict)
            X.append(feat_dict)
        if Y_col:
            Y = np.array(table[Y_col].astype(str))
            return X, Y
        else:
            return X, None

    def _extract_protein_features(self, table, Y_col=None):
        """
        Extracts features on protein level
        given
        table:  a table with pfam domains and class lables
        Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
        """
        X = []
        Y = []
        for prot, tbl in table.groupby(prot_col, sort=False):
            feat_dict = dict()
            for _, row in tbl.iterrows():
                feat_dict = self._make_feature_dict(row, feat_dict)
            X.append(feat_dict)
            if Y_col:
                Y.append(str(tbl[Y_col].values[0]))
        if Y_col:
            return X, Y
        else:
            return X, None

    def _extract_overlapping_features(self, table, Y_col=None):
        """
        Prepares class labels Y and features from a table
        given
        table:  a table with pfam domains and class lables
        Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
        """
        X = []
        for idx, _ in table.iterrows():
            wind = table.iloc[idx - overlap : idx + overlap + 1]
            feat_dict = dict()
            for _, row in wind.iterrows():
                feat_dict = self._make_feature_dict(row, feat_dict)
            X.append(feat_dict)
        if Y_col:
            Y = np.array(table[Y_col].astype(str))
            return X, Y
        else:
            return X, None

    def _make_feature_dict(self, row, feat_dict=dict()):
        """
        Constructs a dict with feature:value pairs from
        row: input row or dict
        self.features: either name of feature or name of column in row/dict
        self.weights: either numerical weight or name of column in row/dict
        """
        for f, w in zip(self.features, self.weights):
            if isinstance(w, numbers.Number):
                feat_dict[row[f]] = w
                continue
            try:
                feat_dict[row[f]] = row[w]
            except KeyError:
                feat_dict[f] = row[w]
        return feat_dict
