"""BGC prediction using a conditional random field.
"""

import csv
import copy
import functools
import hashlib
import itertools
import io
import math
import numbers
import operator
import os
import pickle
import random
import textwrap
import typing
import warnings
from typing import (
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
)

import numpy
import tqdm
import sklearn_crfsuite
import sklearn.model_selection
import sklearn.preprocessing

from .._meta import sliding_window
from ..model import Gene
from . import features
from .cv import LeaveOneGroupOut
from .select import fisher_significance

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

__all__ = ["ClusterCRF"]


class ClusterCRF(object):
    """A wrapper for `sklearn_crfsuite.CRF` to work with the GECCO data model.
    """

    @classmethod
    def trained(cls, model_path: Optional[str] = None) -> "ClusterCRF":
        """Create a new pre-trained `ClusterCRF` instance from a model path.

        Arguments:
            model_path (`str`, optional): The path to the model directory
                obtained with the ``gecco train`` command. If `None` given,
                use the embedded model.

        Returns:
            `~gecco.crf.ClusterCRF`: A CRF model that can be used to perform
            predictions without training first.

        Raises:
            `ValueError`: If the model data does not match its hash.

        """
        # get the path to the pickled model and read its signature file
        if model_path is not None:
            pkl_file = open(os.path.join(model_path, "model.pkl"), "rb")
            md5_file = open(os.path.join(model_path, "model.pkl.md5"))
        else:
            pkl_file = importlib_resources.open_binary(__name__, "model.pkl")
            md5_file = importlib_resources.open_text(__name__, "model.pkl.md5")
        with md5_file as sig:
            signature = sig.read().strip()

        # check the file content matches its MD5 hashsum
        hasher = hashlib.md5()
        with pkl_file as bin:
            read = functools.partial(bin.read, io.DEFAULT_BUFFER_SIZE)
            for chunk in iter(read, b""):
                hasher.update(typing.cast(bytes, chunk))
            if hasher.hexdigest().upper() != signature.upper():
                raise ValueError("MD5 hash of model data does not match signature")
            pkl_file.seek(0)
            return pickle.load(bin)  # type: ignore

    def __init__(
        self,
        feature_type: str = "protein",
        algorithm: str = "lbfgs",
        window_size: int = 5,
        window_step: int = 1,
        **kwargs: Dict[str, object],
    ) -> None:
        """Create a new `ClusterCRF` instance.

        Arguments:
            feature_type (`str`): Defines how features should be extracted.
                Should be either *domain* or *protein*.
            algorithm (`str`): The optimization algorithm for the model. See
                https://sklearn-crfsuite.readthedocs.io/en/latest/api.html
                for available values.
            window_size (`int`):

        Any additional keyword argument is passed as-is to the internal
        `~sklearn_crfsuite.CRF` constructor.

        Raises:
            `ValueError`: if ``feature_type`` has an invalid value.
            `TypeError`: if one of the ``*_columns`` argument is not iterable.

        """
        if feature_type not in {"protein", "domain"}:
            raise ValueError(f"invalid feature type: {feature_type!r}")
        if window_size <= 0:
            raise ValueError("Window size must be strictly positive")
        if window_step <= 0 or window_step > window_size:
            raise ValueError("Window step must be strictly positive and under `window_size`")

        self.feature_type = feature_type
        self.window_size = window_size
        self.window_step = window_step
        self.algorithm = algorithm
        self.significance: Optional[Dict[str, float]] = None
        self.significant_features: Optional[FrozenSet[str]] = None
        self.model = sklearn_crfsuite.CRF(
            algorithm=algorithm,
            # all_possible_transitions=True,
            # all_possible_states=True,
            **kwargs,
        )

    def predict_probabilities(self, genes: Iterable[Gene], *, cpus: Optional[int] = None) -> List[Gene]:
        """Predict how likely each given gene is part of a gene cluster.
        """
        # select the feature extraction method
        if self.feature_type == "protein":
            extract_features = features.extract_features_protein
            annotate_probabilities = features.annotate_probabilities_protein
        elif self.feature_type == "domain":
            extract_features = features.extract_features_domain
            annotate_probabilities = features.annotate_probabilities_domain
        else:
            raise ValueError(f"invalid feature type: {self.feature_type!r}")

        # sort genes by sequence id and domains inside genes by coordinate
        genes = sorted(genes, key=operator.attrgetter("source.id"))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start"))

        # predict sequence by sequence
        predicted = []
        for _, group in itertools.groupby(genes, key=operator.attrgetter("source.id")):
            # extract features
            sequence: List[Gene] = sorted(group, key=operator.attrgetter("start"))
            feats: List[Dict[str, bool]] = extract_features(sequence)
            # ignore sequences too small with a warning
            if len(feats) < self.window_size:
                warnings.warn(
                    f"Contig {sequence[0].source.id!r} does not contain enough"
                    f" genes for sliding window of size {self.window_size}"
                )
                continue
            # predict marginals over a sliding window, storing maximum probabilities
            probabilities = numpy.zeros(len(sequence))
            for win in sliding_window(len(feats), self.window_size, self.window_step):
                marginals = [p['1'] for p in self.model.predict_marginals_single(feats[win])]
                numpy.maximum(probabilities[win], marginals, out=probabilities[win])
            # label genes with maximal probabilities
            predicted.extend(annotate_probabilities(sequence, probabilities))

        # return the genes that were passed as input but now having BGC
        return predicted

    def fit(
        self,
        genes: Iterable[Gene],
        select: Optional[float] = None,
        shuffle: bool = True,
        *,
        cpus: Optional[int] = None,
        correction_method: Optional[str] = None,
    ) -> None:
        _cpus = os.cpu_count() if not cpus else cpus
        # select the feature extraction method
        if self.feature_type == "protein":
            extract_features = features.extract_features_protein
            extract_labels = features.extract_labels_protein
        elif self.feature_type == "domain":
            extract_features = features.extract_features_domain
            extract_labels = features.extract_labels_domain
        else:
            raise ValueError(f"invalid feature type: {self.feature_type!r}")

        # sort genes by sequence id and domains inside genes by coordinate
        genes = sorted(genes, key=operator.attrgetter("source.id"))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start"))

        # perform feature selection
        if select is not None:
            if select <= 0 or select > 1:
                raise ValueError(f"invalid value for select: {select}")
            # find most significant features
            self.significance = sig = fisher_significance(
                (gene.protein for gene in genes),
                correction_method=correction_method,
            )
            sorted_sig = sorted(sig, key=sig.get)[:int(select*len(sig))]
            self.significant_features = frozenset(sorted_sig)
            # check that we don't select "random" domains
            if sig[sorted_sig[-1]] == 1.0:
                warnings.warn(
                    "Selected features still include domains with a p-value "
                    "of 1, consider reducing the selected fraction.",
                    UserWarning
                )
            # remove non significant domains
            genes = [
                gene.with_protein(
                    gene.protein.with_domains([
                        domain for domain in gene.protein.domains
                        if domain.name in self.significant_features
                    ])
                )
                for gene in genes
            ]

        # group input genes by sequence
        groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
        sequences = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]
        # shuffle sequences if requested
        if shuffle:
            random.shuffle(sequences)

        # build training instances
        training_features, training_labels = [], []
        for sequence in sequences:
            # extract features and labels
            feats: List[Dict[str, bool]] = extract_features(sequence)
            labels: List[str] = extract_labels(sequence)
            # check we have as many observations as we have labels
            if len(feats) != len(labels):
                raise ValueError("different number of features and labels found, something is wrong")
            # check we have enough genes for the desired sliding window
            if len(feats) < self.window_size:
                raise ValueError(f"{sequence[0].source.id!r} has not enough observations ({len(feats)}) for requested window size ({self.window_size})")
            # record every window in the sequence
            for win in sliding_window(len(feats), self.window_size, self.window_step):
                training_features.append(feats[win])
                training_labels.append(labels[win])

        # check labels
        if all(label == "1" for y in training_labels for label in y):
            raise ValueError("only positives labels found, something is wrong.")
        elif all(label == "0" for y in training_labels for label in y):
            raise ValueError("only negative labels found, something is wrong.")

        # fit the model
        self.model.fit(training_features, training_labels)
