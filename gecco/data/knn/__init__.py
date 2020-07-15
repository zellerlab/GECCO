import csv
import itertools
import gzip
import pkg_resources
import typing

import numpy


class TrainingMatrix(typing.NamedTuple):
    compositions: numpy.ndarray
    labels: numpy.ndarray
    types: numpy.ndarray
    domains: numpy.ndarray


def load_training_matrix() -> TrainingMatrix:
    with pkg_resources.resource_stream(__name__, "training_matrix.tsv.gz") as bin:
        stream = gzip.open(bin, mode="rt")
        reader = csv.reader(stream, dialect="excel-tab")
        compositions, labels, types, domains = [], [], [], next(reader)[2:]
        for line in reader:
            labels.append(line[0])
            types.append(line[1])
            compositions.append(numpy.array(line[2:]).astype(float))
        return TrainingMatrix(
            compositions=numpy.array(compositions),
            labels=numpy.array(labels),
            types=numpy.array(types),
            domains=numpy.array(domains)
        )
