"""Supervised classifier to predict the biosynthetic type of a cluster.
"""

import csv
import functools
import operator
import os
import typing
from typing import Callable, Dict, List, Iterable, Optional, Sequence, Tuple, Union

import numpy
import pkg_resources
import scipy.sparse
import sklearn.ensemble
import sklearn.preprocessing

from ..model import ProductType, Cluster


class TypeBinarizer(sklearn.preprocessing.MultiLabelBinarizer):

    def __init__(self):
        self.classes_ = sorted(x for x in ProductType.__members__.values() if x)
        super().__init__(self.classes_)

    def transform(self, y: Iterable[ProductType]) -> Iterable[Iterable[int]]:
        matrix = numpy.zeros((len(y), len(self.classes)))
        for i, label in enumerate(y):
            for j, cls in enumerate(self.classes_):
                matrix[i, j] = int(label & cls)
        return matrix

    def inverse_transform(self, yt: Iterable[Iterable[int]]) -> Iterable[ProductType]:
        classes = []
        for y in yt:
            filtered = [cls for i, cls in enumerate(self.classes_) if y[i]]
            classes.append(functools.reduce(operator.or_, filtered))
        return classes



class TypeClassifier(object):

    @classmethod
    def trained(cls, model_path: Optional[str] = None) -> "TypeClassifier":
        """Create a new `TypeClassifier` pre-trained with embedded data.

        Arguments:
            model_path (`str`, optional): The path to the model directory
                obtained with the ``gecco train`` command. If `None` given,
                use the embedded training data.

        Returns:
            `~gecco.types.TypeClassifier`: A random forest model that can be
            used to perform BGC type predictions without training first.

        """

        if model_path is not None:
            doms_path = os.path.join(model_path, "domains.tsv")
            typs_path = os.path.join(model_path, "types.tsv")
            comp_path = os.path.join(model_path, "compositions.npz")
        else:
            doms_path = pkg_resources.resource_filename(__name__, "domains.tsv")
            typs_path = pkg_resources.resource_filename(__name__, "types.tsv")
            comp_path = pkg_resources.resource_filename(__name__, "compositions.npz")

        compositions = scipy.sparse.load_npz(comp_path)
        with open(doms_path, "r") as src:
            domains = [ line.strip() for line in src ]
        with open(typs_path, "r") as src:
            types = [
                [ProductType.__members__[ty] for ty in tys.split(";")]
                for tys in (line.split("\t")[1].strip() for line in src)
            ]
            compact_types = [
                functools.reduce(operator.or_, tys) for tys in types
            ]

        classifier = cls(random_state=0)
        types_bin = classifier.binarizer.transform(compact_types)
        classifier.model.fit(compositions, y=types_bin)
        classifier.model.attributes_ = domains
        return classifier

    def __init__(self, **kwargs: object) -> None:
        """Instantiate a new type classifier.

        Keyword Arguments:
            Any additional keyword argument is passed as argument to the
            internal `~sklearn.ensemble.RandomForestClassifier` constructor.

        """
        self.model = sklearn.ensemble.RandomForestClassifier(**kwargs)
        self.binarizer = TypeBinarizer()

    _S = typing.TypeVar("_S", bound=Sequence["Cluster"])

    def predict_types(self, clusters: "_S") -> "_S":
        """Predict types for each of the given clusters.
        """
        # extract domain compositions from input clusters
        comps = numpy.array([c.domain_composition(self.model.attributes_) for c in clusters])

        # only take 'positive' probabilities for each class
        probas = numpy.array(self.model.predict_proba(comps))[:, :, 1].transpose()
        preds = self.binarizer.inverse_transform(probas > 0.5)

        # annotate the input clusters
        results = zip(typing.cast(Iterable["Cluster"], clusters), probas, preds)
        for cluster, proba, ty in results:
            cluster.type = ty
            cluster.type_probabilities = dict(zip(self.binarizer.classes_, proba))
        return clusters
