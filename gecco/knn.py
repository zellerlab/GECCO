import functools
import typing
from typing import Callable, Dict, Optional, Union

import sklearn.neighbors
import scipy.spatial.distance

if typing.TYPE_CHECKING:
    import numpy
    _Metric = typing.Callable[["numpy.ndarray", "numpy.ndarray"], "numpy.ndarray"]


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

    @classmethod
    def _get_metric(cls, name: str) -> "_Metric":
        if name == "jensenshannon":
            return scipy.spatial.distance.jensenshannon  # type: ignore
        elif name == "tanimoto":
            # NB(@althonos): Tanimoto distance seems to be mostly for boolean
            #                vectors, not probability vectors.
            return lambda p,q: p*q / (p - q)**2  # type: ignore
        else:
            raise ValueError(f"unexpected metric: {name!r}")

    def __init__(
            self,
            metric: Union[str, "_Metric"] = "jensenshannon",
            **kwargs: object
    ) -> None:
        """Create a new classifier.

        Arguments:
            metric (`str` or `function`): The distance metric to use with the
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
            self.metric = self._get_metric(metric)
        else:
            self.metric = metric
        self.knn = sklearn.neighbors.KNeighborsClassifier(
            metric=self.metric,
            **kwargs
        )

    def fit_predict(self, train_matrix, new_matrix, y):
        """Fit the model and immediately produce a prediction.
        """
        self.knn.fit(train_matrix, y=y)
        pred = self.knn.predict(new_matrix)
        proba = [max(p) for p in self.knn.predict_proba(new_matrix)]
        return list(zip(pred, proba))
