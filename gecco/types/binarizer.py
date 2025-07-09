import typing
from typing import List, Iterable

import numpy
import sklearn.preprocessing

from ..model import ClusterType

if typing.TYPE_CHECKING:
    from numpy.typing import NDArray


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
