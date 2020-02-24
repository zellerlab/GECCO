import typing

import numpy
from scipy.spatial.distance import jensenshannon
from sklearn.neighbors import KNeighborsClassifier


_Metric = typing.Callable[[numpy.ndarray, numpy.ndarray], numpy.ndarray]


class ClusterKNN(object):
    """A kNN classifier to predict cluster types based on domain similarity.

    Essentially a wrapper around `sklearn.neighbors.KNeighborsClassifier`
    provided for convenience.

    Attributes:
        metric (str): The name of the distance metric used by the classifier.
        dist (function): The actual function performing the distance
            computation.
        knn (`sklearn.neighbors.KNeighborsClassifier`): The internal kNN
            classifier used for the prediction
    """

    _METRICS: typing.Dict[str, _Metric] = {
        # NB(@althonos): Tanimoto distance seems to be mostly for boolean
        #                vectors, not probability vectors.
        "tanimoto": lambda p,q: p*q / (p - q)**2,
        "jensenshannon": jensenshannon,
    }

    def __init__(
            self,
            metric: typing.Union[str, _Metric] = "jensenshannon",
            **kwargs: object
    ):
        """Create a new classifier.

        Arguments:
            metric (str or function): The distance metric to use with the
                classifier. Either given a metric name (such as
                ``jensenshannon``, the default) or a callable that takes
                two vectors.

        Keyword Arguments:
            Any additional keyword argument is passed as argument to the
            internal `~sklearn.neighbors.KNeighborsClassifier` constructor.

        Raises:
            ValueError: when an unknown distance metric is given.
        """

        if isinstance(metric, str):
            self.metric = self._METRICS.get(metric)
        else:
            self.metric = metric
        if self.metric is None:
            raise ValueError(f"unexpected metric: {metric!r}")
        self.knn = KNeighborsClassifier(metric=self.metric, **kwargs)

    def fit_predict(self, train_matrix, new_matrix, y):
        """Fit the model and immediately produce a prediction.
        """
        self.knn.fit(train_matrix, y=y)
        pred = self.knn.predict(new_matrix)
        proba = [max(p) for p in self.knn.predict_proba(new_matrix)]
        return list(zip(pred, proba))
