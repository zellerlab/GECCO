import numpy as np
import pandas as pd
from sklearn.manifold import TSNE, MDS
from sklearn.neighbors import KNeighborsClassifier
from scipy.spatial.distance import jensenshannon


class ClusterKNN(object):
    """
    Predicts cluster types based on a kNN Classifier and plots results.
    Essentially a wrapper around sklearns KNeighborsClassifier
    (MDS and TSNE maybe later)
    """

    _METRICS = {
        # Doesn't work, really
        # NB(@althonos): Tanimoto distance seems to be mostly for boolean
        #                vectors, not probability vectors.
        "tanimoto": lambda p,q: p*q / (p - q)**2,
        "jensenshannon": jensenshannon,
    }

    def __init__(self, metric: str = "jensenshannon", **kwargs):
        self.metric = metric
        self.dist = dist = self._METRICS.get(metric)
        if dist is None:
            raise ValueError(f"unexpected metric: {metric!r}")
        self.knn = KNeighborsClassifier(metric=dist, **kwargs)

    def fit_predict(self, train_matrix, new_matrix, y):
        self.knn.fit(train_matrix, y=y)
        pred = self.knn.predict(new_matrix)
        proba = [max(p) for p in self.knn.predict_proba(new_matrix)]
        return list(zip(pred, proba))
