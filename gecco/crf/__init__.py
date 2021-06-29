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
import sklearn_crfsuite
import sklearn.model_selection
import sklearn.preprocessing

from ..model import Gene
from .._meta import OrderedPoolWrapper
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
        self.significance: Optional[Dict[str, float]] = None
        self.significant_features: Optional[FrozenSet[str]] = None
        self.model = sklearn_crfsuite.CRF(
            algorithm=algorithm,
            all_possible_transitions=True,
            all_possible_states=True,
            **kwargs,
        )

    def predict_probabilities(self, genes: Iterable[Gene], *, cpus: Optional[int] = None) -> List[Gene]:
        """Predict how likely each given gene is part of a gene cluster.
        """
        _cpus = os.cpu_count() if not cpus else cpus
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
        with OrderedPoolWrapper(self.pool_factory(_cpus)) as pool:
            # extract features in parallel and predict cluster probabilities
            marginals = self.model.predict_marginals(pool.map(extract, seqs))
            # Annotate the genes with the predicted probabilities
            annotated_seqs = pool.starmap(annotate, zip(seqs, marginals))

        # return the genes that were passed as input but now having BGC
        # probabilities set
        return list(itertools.chain.from_iterable(annotated_seqs)) # type: ignore

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
        if self.feature_type == "group":
            extract_features = features.extract_features_group
            extract_labels = features.extract_labels_group
        elif self.feature_type == "single":
            extract_features = features.extract_features_single
            extract_labels = features.extract_labels_single
        elif self.feature_type == "overlap":
            raise NotImplementedError("todo: `features.extract_features_overlap`")
        else:
            raise ValueError("invalid feature type")

        # group input genes by sequence
        groups = itertools.groupby(genes, key=operator.attrgetter("source.id"))
        seqs = [sorted(group, key=operator.attrgetter("start")) for _, group in groups]

        # shuffle sequences
        if shuffle:
            random.shuffle(seqs)

        # perform feature selection
        if select is not None:
            if select <= 0 or select > 1:
                raise ValueError(f"invalid value for select: {select}")
            # find most significant features
            self.significance = sig = fisher_significance(
                (g.protein for seq in seqs for g in seq),
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
            for i, seq in enumerate(seqs):
                for j, gene in enumerate(seq):
                    seqs[i][j] = gene = copy.deepcopy(gene)
                    gene.protein.domains = [
                        domain for domain in gene.protein.domains
                        if domain.name in self.significant_features
                    ]

        # proces each sequence / group in parallel
        with OrderedPoolWrapper(self.pool_factory(_cpus)) as pool:
            # extract features in parallel and predict cluster probabilities
            X = pool.map(extract_features, seqs)
            Y = pool.map(extract_labels, seqs)

        # check labels
        if all(y == "1" for y in Y):
            raise ValueError("only positives labels found, something is wrong.")
        elif all(y == "0" for y in Y):
            raise ValueError("only negative labels found, something is wrong.")

        # fit the model
        self.model.fit(X, Y)
