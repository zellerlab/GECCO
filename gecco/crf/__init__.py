"""BGC prediction using a conditional random field.
"""

import csv
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
from multiprocessing.pool import Pool
from typing import (
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Tuple,
    Type,
    Union,
)

import numpy
import tqdm
import pkg_resources
import sklearn_crfsuite
import sklearn.model_selection
import sklearn.preprocessing

from ..model import Gene
from .._meta import OrderedPoolWrapper
from . import features
from .cv import LeaveOneGroupOut
from .select import fisher_significance


class ClusterCRF(object):

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
            pkl_path = os.path.join(model_path, "model.pkl")
            md5_path = os.path.join(model_path, "model.pkl.md5")
        else:
            pkl_path = pkg_resources.resource_filename(__name__, "model.pkl")
            md5_path = pkg_resources.resource_filename(__name__, "model.pkl.md5")
        with open(md5_path) as sig:
            signature = sig.read().strip()

        # check the file content matches its MD5 hashsum
        hasher = hashlib.md5()
        with open(pkl_path, "rb") as bin:
            read = functools.partial(bin.read, io.DEFAULT_BUFFER_SIZE)
            for chunk in iter(read, b""):
                hasher.update(typing.cast(bytes, chunk))
        if hasher.hexdigest().upper() != signature.upper():
            raise ValueError("MD5 hash of model data does not match signature")

        # load the pickled model if the data matches the hashsum
        with open(pkl_path, "rb") as bin:
            return pickle.load(bin)  # type: ignore

    def __init__(
        self,
        feature_type: str = "protein",
        algorithm: str = "lbfgs",
        overlap: int = 2,
        pool_factory: Union[Type[Pool], Callable[[Optional[int]], Pool]] = Pool,
        **kwargs: Dict[str, object],
    ) -> None:
        """Create a new `ClusterCRF` instance.

        Arguments:
            feature_type (`str`): Defines how features should be extracted.
                Should be either *domain*, *protein*, or *overlap*.
            algorithm (`str`): The optimization algorithm for the model. See
                https://sklearn-crfsuite.readthedocs.io/en/latest/api.html
                for available values.
            overlap (`int`): In case of ``feature_type = "overlap"``, defines
                the sliding window size to use. The resulting window width is
                ``2*overlap+1``.
            pool_factory (`multiprocessing.pool.Pool` subclass, or callable):
                The factory to use to create a new pool instance for methods
                that can perform operations in parallel. *It is called with
                a single argument which is the number of workers to create,
                or* `None` *to create a much workers as there are CPUs.*

        Any additional keyword argument is passed as-is to the internal
        `~sklearn_crfsuite.CRF` constructor.

        Raises:
            `ValueError`: if ``feature_type`` has an invalid value.
            `TypeError`: if one of the ``*_columns`` argument is not iterable.

        """
        if feature_type not in {"single", "overlap", "group"}:
            raise ValueError(f"unexpected feature type: {feature_type!r}")

        self.feature_type: str = feature_type
        self.overlap: int = overlap
        self.algorithm = algorithm
        self.pool_factory = pool_factory
        self.significant_features: Dict[str, FrozenSet[str]] = {}
        self.model = sklearn_crfsuite.CRF(
            algorithm=algorithm,
            all_possible_transitions=True,
            all_possible_states=True,
            **kwargs,
        )

    def predict_probabilities(self, genes: Iterable[Gene], *, jobs: Optional[int] = None) -> List[Gene]:
        """Predict
        """
        # group input genes by sequence
        groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
        seqs = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]

        # select the feature extraction method
        if self.feature_type == "group":
            extract = features.extract_features_group
            annotate = features.annotate_probabilities_group
        elif self.feature_type == "single":
            extract = features.extract_features_single
            annotate = features.annotate_probabilities_single
        elif self.feature_type == "overlap":
            raise NotImplementedError("todo: `features.extract_features_overlap`")
        else:
            raise ValueError("invalid feature type")

        # proces each sequence / group in parallel
        with OrderedPoolWrapper(self.pool_factory(jobs)) as pool:
            # extract features in parallel and predict cluster probabilities
            marginals = self.model.predict_marginals(pool.map(extract, seqs))
            # Annotate the genes with the predicted probabilities
            annotated_genes = pool.starmap(annotate, zip(seqs, marginals))

        # return the genes that were passed as input but now having BGC
        # probabilities set
        return list(itertools.chain.from_iterable(annotated_genes))
