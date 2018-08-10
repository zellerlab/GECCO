import numpy as np
import pandas as pd
from sklearn.manifold import TSNE, MDS
from sklearn.neighbors import KNeighborsClassifier
from orion.utils import jsd_pairwise, tanimoto_pairwise

class ClusterKNN(object):
    """
    Predicts cluster types based on a kNN Classifier and plots results.
    Essentially a wrapper around sklearns KNeighborsClassifier
    (MDS and TSNE maybe later)
    """

    def __init__(self, metric="jsd", **kwargs):
        self.metric = metric

        if self.metric == "jsd":
            self.dist = jsd_pairwise

        if self.metric == "tanimoto":
            self.dist = tanimoto_pairwise

        self.knn = KNeighborsClassifier(
            metric = self.dist,
            **kwargs
        )

    def fit_predict(self, train_matrix, new_matrix, y):
        self.knn.fit(train_matrix, y=y)
        pred = self.knn.predict(new_matrix)
        proba = [max(p) for p in self.knn.predict_proba(new_matrix)]
        return list(zip(pred, proba))
