import csv
import io
import itertools
import gzip
import pkg_resources
import typing

import numpy
import scipy.sparse


class TrainingMatrix(typing.NamedTuple):
    compositions: numpy.ndarray
    labels: numpy.ndarray
    types: numpy.ndarray
    domains: numpy.ndarray


def load_training_matrix() -> TrainingMatrix:
    with pkg_resources.resource_stream(__name__, "domains.tsv") as f:
        domains = f.read().decode("ascii").splitlines()
    with pkg_resources.resource_stream(__name__, "types.tsv") as f:
        reader = csv.reader(io.TextIOWrapper(f), dialect="excel-tab")
        labels, types = [], []
        for row in reader:
            labels.append(row[0])
            types.append(row[1])
    with pkg_resources.resource_stream(__name__, "compositions.npz") as f:
        compositions = scipy.sparse.load_npz(f)
    return TrainingMatrix(
        compositions=compositions.toarray(),
        labels=numpy.array(labels),
        types=numpy.array(types),
        domains=numpy.array(domains)
    )
