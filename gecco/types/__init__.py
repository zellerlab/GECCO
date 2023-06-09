"""Supervised classifier to predict the type of a cluster.
"""

import csv
import functools
import operator
import os
import typing
import warnings
from typing import (
    BinaryIO,
    Callable,
    ContextManager,
    Dict,
    List,
    Iterable,
    Optional,
    Sequence,
    TextIO,
    Tuple,
    Union
)

import numpy
import scipy.sparse
import sklearn.ensemble
import sklearn.preprocessing

from ..model import ClusterType, Cluster

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files  # type: ignore

if typing.TYPE_CHECKING:
    from numpy.typing import NDArray


__all__ = ["TypeBinarizer", "TypeClassifier"]


class TypeBinarizer(sklearn.preprocessing.MultiLabelBinarizer):
    """A `MultiLabelBinarizer` working with `ClusterType` instances.
    """

    def __init__(self, classes: List[str], **kwargs: object):
        self.classes_ = classes
        super().__init__(classes=classes, **kwargs)

    def transform(self, y: List[ClusterType]) -> Iterable[Iterable[int]]:
        matrix = numpy.zeros((len(y), len(self.classes_)))
        for i, label in enumerate(y):
            for j, cls in enumerate(self.classes_):
                matrix[i, j] = cls in label.names
        return matrix

    def inverse_transform(self, yt: "NDArray[numpy.bool_]") -> Iterable[ClusterType]:
        classes = []
        for y in yt:
            filtered = (cls for i, cls in enumerate(self.classes_) if y[i])
            classes.append(ClusterType(*filtered))
        return classes


class TypeClassifier(object):
    """A wrapper to predict the type of a `~gecco.model.Cluster`.
    """

    @classmethod
    def trained(cls, model_path: Optional[str] = None) -> "TypeClassifier":
        """Create a new `TypeClassifier` pre-trained with embedded data.

        Arguments:
            model_path (`str`, optional): The path to the model directory
                obtained with the ``gecco train`` command. If `None` given,
                use the embedded training data.

        Returns:
            `~gecco.types.TypeClassifier`: A random forest model that can be
            used to perform cluster type predictions without training first.

        """

        if model_path is not None:
            doms_file: ContextManager[TextIO] = open(os.path.join(model_path, "domains.tsv"))
            typs_file: ContextManager[TextIO] = open(os.path.join(model_path, "types.tsv"))
            comp_file: ContextManager[BinaryIO] = open(os.path.join(model_path, "compositions.npz"), "rb")
        else:
            doms_file = files(__name__).joinpath("domains.tsv").open()
            typs_file = files(__name__).joinpath("types.tsv").open()
            comp_file = files(__name__).joinpath("compositions.npz").open("rb")

        with comp_file as comp_src:
            compositions = scipy.sparse.load_npz(comp_src)
        with doms_file as doms_src:
            domains = [ line.strip() for line in doms_src ]
        with typs_file as typs_src:
            types = []
            unique_types = set()
            for line in typs_src:
                unpacked = set()
                for ty in filter(None, line.split("\t")[1].strip().split(";")):
                    unpacked.add(ty)
                    unique_types.add(ty)
                types.append(ClusterType(*unpacked))

        classifier = cls(classes=sorted(unique_types), random_state=0)
        types_bin = classifier.binarizer.transform(types)
        if len(classifier.classes_) > 1:
            classifier.model.fit(compositions, y=types_bin)
        classifier.model.attributes_ = domains
        return classifier

    def __init__(self, classes: Iterable[str] = (), **kwargs: object) -> None:
        """Instantiate a new type classifier.

        Keyword Arguments:
            Any additional keyword argument is passed as argument to the
            internal `~sklearn.ensemble.RandomForestClassifier` constructor.

        """
        self.model = sklearn.ensemble.RandomForestClassifier(**kwargs)
        self.binarizer = TypeBinarizer(list(classes))

    @property
    def classes_(self) -> List[str]:
        return self.binarizer.classes_

    _S = typing.TypeVar("_S", bound=Sequence["Cluster"])

    def predict_types(self, clusters: "_S") -> "_S":
        """Predict types for each of the given clusters.
        """
        # extract domain compositions from input clusters
        comps = numpy.array([c.domain_composition(self.model.attributes_) for c in clusters])

        # predict type probabilities with the internal classifier
        probas = self.model.predict_proba(comps)

        # extract only the *positive* probabilites and translate them to proper
        # type predictions using the binarizer
        if len(comps) == 1:
            posit = numpy.array([[1 - cls[0][0] for cls in probas]])
        else:
            posit = 1 - numpy.array(probas)[:, :, 0].transpose()

        # translate probabilities into product type predictions
        types = self.binarizer.inverse_transform(posit > 0.5)

        # annotate the input clusters
        results = zip(typing.cast(Iterable["Cluster"], clusters), posit, types)
        for cluster, proba, ty in results:
            cluster.type = ty
            cluster.type_probabilities = dict(zip(self.binarizer.classes_, proba))
        return clusters
