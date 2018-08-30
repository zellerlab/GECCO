import math
import numbers
import numpy as np
import pandas as pd
from math import ceil
from itertools import zip_longest
from collections import Iterable
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.base import TransformerMixin, BaseEstimator

# CLASS
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


# FUNC
def flatten(l):
    """Flattens list of arbitraty depth to generator."""
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def truncate(df, length, Y_col="BGC", grouping="protein_id"):

    df0 = df[df[Y_col] == 0]
    df1 = df[df[Y_col] == 1]
    df0 = [df for _, df in df0.groupby(grouping, sort=False)]
    trunc_len = int((len(df0) - 2 * length) / 2)

    try:
        df_trunc = pd.concat(df0[trunc_len : -trunc_len])
    except ValueError as err:
        df_trunc = pd.concat(df0)

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

def prepare_sample(table, X_col, Y_col):
    """
    Prepares features X (Pfam domains) and class labels Y from a table
    given
    table:  a table with pfam domains and class lables
    X_col:  the column containing the Pfam domain annotation
    Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
    """
    X = table[X_col]
    X = np.array(X).reshape(-1, 1)
    Y = np.array(table[Y_col])
    return X, Y

def prepare_split_samples(table_list, split_col, X_col="pfam",
        Y_col="BGC", ascending=True):
    """
    Splits tables into separate samples
    given
    table_list:  a list of paths to sample tabels
    split_col: column to split into separate samples
    separator:  the separator to read these tables
    X_col:  the column containing the Pfam domain annotation
    Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
    sort_col: the column to sort the values
    ascending: sort ascending vs. descending.
    """
    X = []
    Y = []
    for tbl in list(table_list):
        samples =  [prepare_sample(s, X_col, Y_col) for _, s in total.groupby(split_col)]
        X += [x for x, _ in samples]
        Y += [y for _, y in samples]
    return np.array(X), np.array(Y)

def filter_repeats(pfam_df):
    repeats = set(['PF07721', 'PF05593', 'PF07719', 'PF00515', 'PF00132', 'PF03130',
        'PF01839', 'PF01816', 'PF07720', 'PF00400', 'PF05594', 'PF07661', 'PF02985',
        'PF06049', 'PF08238', 'PF06696', 'PF00353', 'PF02412', 'PF00023', 'PF02071',
        'PF03991', 'PF01469', 'PF07676', 'PF00514', 'PF00904', 'PF07634', 'PF02370',
        'PF03335', 'PF01851', 'PF04728', 'PF06715', 'PF03373', 'PF04680', 'PF00805',
        'PF04508', 'PF07918', 'PF01535', 'PF01011', 'PF05017', 'PF06671', 'PF00818',
        'PF03406', 'PF00399', 'PF09373', 'PF01744', 'PF01436', 'PF01239', 'PF05906',
        'PF03729', 'PF00404', 'PF04022', 'PF02363', 'PF02524', 'PF07981', 'PF02095',
        'PF00414', 'PF00560', 'PF05001', 'PF02162', 'PF01473', 'PF05465', 'PF02493',
        'PF03578', 'PF08043', 'PF06392', 'PF07142', 'PF08309', 'PF02184'])
    out_df = pfam_df[~pfam_df["pfam"].isin(repeats)]
    return out_df

def extract_features(table, Y_col=None, feature_col=[], weight_col=[]):
    """
    Prepares class labels Y and features from a table
    given
    table:  a table with pfam domains and class lables
    Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
    feature_col: either name of a column with a categorical feature or the name of a numerical feature
    weight_col: either name of a column with numerical values or a numerical weight
    """
    X = []
    for _, row in table.iterrows():
        feat_dict = dict()
        feat_dict = make_feature_dict(row, feature_col, weight_col, feat_dict)
        X.append(feat_dict)
    if Y_col:
        Y = np.array(table[Y_col].astype(str))
        return X, Y
    else:
        return X, None

def extract_protein_features(table, Y_col=None, feature_col="pfam", prot_col="protein_id",
        weight_col="rev_i_Evalue"):
    """
    Extracts features on protein level
    given
    table:  a table with pfam domains and class lables
    Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
    feature_col: either name of a column with a categorical feature or the name of a numerical feature
    weight_col: either name of a column with numerical values or a numerical weight
    """
    X = []
    Y = []
    for prot, tbl in table.groupby(prot_col, sort=False):
        feat_dict = dict()
        for _, row in tbl.iterrows():
            feat_dict = make_feature_dict(row, feature_col, weight_col, feat_dict)
        X.append(feat_dict)
        if Y_col:
            Y.append(str(tbl[Y_col].values[0]))
    if Y_col:
        return X, Y
    else:
        return X, None

def extract_overlapping_features(table, Y_col=None, feature_col=[], weight_col=[],
        overlap=1):
    """
    Prepares class labels Y and features from a table
    given
    table:  a table with pfam domains and class lables
    Y_col:  the column containing the class labels (e.g. 1 (BGC) and 0 (non-BGC))
    feature_col: either name of a column with a categorical feature or the name of a numerical feature
    weight_col: either name of a column with numerical values or a numerical weight
    """
    X = []
    for idx, _ in table.iterrows():
        wind = table.iloc[idx - overlap : idx + overlap + 1]
        feat_dict = dict()
        for _, row in wind.iterrows():
            feat_dict = make_feature_dict(row, feature_col, weight_col, feat_dict)
        X.append(feat_dict)
    if Y_col:
        Y = np.array(table[Y_col].astype(str))
        return X, Y
    else:
        return X, None

def make_feature_dict(row, features=[], weights=[], feat_dict=dict()):
    """
    Constructs a dict with feature:value pairs from
    row: input row or dict
    features: either name of feature or name of column in row/dict
    weights: either numerical weight or name of column in row/dict
    """
    for f, w in zip(features, weights):
        if isinstance(w, numbers.Number):
            feat_dict[row[f]] = w
            continue
        try:
            feat_dict[row[f]] = row[w]
        except KeyError:
            feat_dict[f] = row[w]
    return feat_dict

def read_to_set(file, split_at="."):
    out_set = set()
    with open(file, "rt") as f:
        for line in f:
            out_set.add(line.strip().split(split_at)[0])

    return out_set

def compute_features(pfam_df, weight_type=None):

    pfam_df = pfam_df.assign(
        pfam = pfam_df["pfam"].str.replace(r"(PF\d+)\.\d+", lambda m: m.group(1)))

    if weight_type == "rev_i_Evalue":
        return pfam_df.assign(rev_i_Evalue = 1 - pfam_df["i_Evalue"])

    if weight_type == "supnorm_Evalue":
        return (pfam_df
            .assign(
                log_Evalue = [min(- np.log10(n), 300) for n in pfam_df["i_Evalue"]]
            )
            .groupby("pfam", sort = False)
            .apply(lambda x: x.assign(
                supnorm_Evalue = (x["log_Evalue"] / x["log_Evalue"].max())
            ))
            .reset_index(drop = True)
        )

    if weight_type == "supnorm_Evalue_prot":
        return (pfam_df
            .assign(
                log_Evalue = [min(- np.log10(n), 300) for n in pfam_df["i_Evalue"]]
            )
            .groupby("pfam", sort = False)
            .apply(lambda x: x.assign(
                supnorm_Evalue = (x["log_Evalue"] / x["log_Evalue"].max())
            ))
            .reset_index(drop = True)
            .groupby("protein_id", sort = False)
            .apply(lambda x: x.assign(
                supnorm_Evalue_prot = x["supnorm_Evalue"] / x["supnorm_Evalue"].sum()
            ))
            .reset_index(drop = True)
            .sort_values("pseudo_pos")
        )

    else:
        return pfam_df
