"""Supervised classifier to predict the biosynthetic class of a putative BGC.
"""

import csv
import functools
import os
import typing
from typing import Callable, Dict, List, Iterable, Optional, Sequence, Tuple, Union

import numpy
import pandas
import pkg_resources
import sklearn.neighbors
import scipy.spatial.distance

from ..model import Cluster

if typing.TYPE_CHECKING:
    import numpy

    _S = typing.TypeVar("_S", bound=Sequence[Cluster])
    _Metric = typing.Callable[["numpy.ndarray", "numpy.ndarray"], "numpy.ndarray"]


class ClusterKNN(object):
    """A kNN classifier to predict cluster types based on domain similarity.

    Essentially a wrapper around `sklearn.neighbors.KNeighborsClassifier`
    provided for convenience.

    Attributes:
        metric (distance): A distance metric used by the classifier. It takes
            two probability vectors as arguments, and returns the distance
            between these two vectors. By default, this is
            `~scipy.spatial.distance.jensenshannon`, but another function such
            as `~scipy.spatial.distance.mahalanobis` can be passed as argument
            to the constructor if needed.
        model (`sklearn.neighbors.KNeighborsClassifier`): The internal kNN
            classifier used for the prediction.

    """

    @classmethod
    def trained(
        cls,
        model_path: Optional[str] = None,
        metric: Union[str, "_Metric"] = "jensenshannon"
    ) -> "ClusterKNN":
        """Create a new `ClusterKNN` instance pre-trained with embedded data.

        Arguments:
            model_path (`str`, optional): The path to the model directory
                obtained with the ``gecco train`` command. If `None` given,
                use the embedded training data.
            metric (`str` or `function`): The distance metric to use with the
                classifier. Either given a metric name (such as
                ``jensenshannon``, the default) or a callable that takes
                two vectors.

        Returns:
            `~gecco.knn.ClusterKNN`: A KNN model that can be used to perform
            predictions without training first.

        """
        if model_path is not None:
            doms_path = os.path.join(model_path, "domains.tsv")
            typs_path = os.path.join(model_path, "types.tsv")
            comp_path = os.path.join(model_path, "compositions.npz")
        else:
            doms_path = pkg_resources.resource_filename(__name__, "domains.tsv")
            typs_path = pkg_resources.resource_filename(__name__, "types.tsv")
            comp_path = pkg_resources.resource_filename(__name__, "compositions.npz")

        compositions = scipy.sparse.load_npz(comp_path).toarray()
        domains = pandas.read_csv(doms_path, header=None, sep="\t")[0].array
        types = pandas.read_csv(typs_path, header=None, sep="\t")[1].array

        knn = cls(metric=metric)
        knn.model.fit(compositions, y=types)
        knn.model.attributes_ = domains
        return knn

    @classmethod
    def _get_metric(cls, name: str) -> "_Metric":
        if name == "jensenshannon":
            return scipy.spatial.distance.jensenshannon  # type: ignore
        elif name == "tanimoto":
            # NB(@althonos): Tanimoto distance seems to be mostly for boolean
            #                vectors, not probability vectors.
            return lambda p, q: p * q / (p - q) ** 2  # type: ignore
        else:
            raise ValueError(f"unknown metric name: {name!r}")

    def __init__(
        self,
        metric: Union[str, "_Metric"] = "jensenshannon",
        **kwargs: object
    ) -> None:
        """Instantiate a new classifier.

        Arguments:
            metric (`str` or `function`): The distance metric to use with the
                classifier. Either given a metric name (such as
                ``jensenshannon``, the default) or a callable that takes
                two vectors.

        Keyword Arguments:
            Any additional keyword argument is passed as argument to the
            internal `~sklearn.neighbors.KNeighborsClassifier` constructor.

        Raises:
            `ValueError`: When an unknown distance metric is given.

        """
        if isinstance(metric, str):
            self.metric = self._get_metric(metric)
        else:
            self.metric = metric
        self.model = sklearn.neighbors.KNeighborsClassifier(metric=self.metric, **kwargs)

    def predict_types(self, clusters: _S) -> _S:
        """Predict types for each of the given clusters.
        """
        comps = [c.domain_composition(self.model.attributes_) for c in clusters]
        probas = self.model.predict_proba(comps)
        for cluster, proba in zip(typing.cast(Iterable[Cluster], clusters), probas):
            imax = numpy.argmax(proba)
            cluster.type = self.model.classes_[imax]
            cluster.type_probability = proba[imax]
        return clusters
