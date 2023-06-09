"""Gene cluster prediction using a conditional random field.
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
    Any,
    BinaryIO,
    Callable,
    ContextManager,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Sequence,
    TextIO,
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
    from importlib.resources import files
except ImportError:
    from importlib_resources import files  # type: ignore

__all__ = ["ClusterCRF"]


class ClusterCRF(object):
    """A wrapper for `sklearn_crfsuite.CRF` to work with the GECCO data model.
    """

    _FILENAME = "model.pkl"

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
            pkl_file: ContextManager[BinaryIO] = open(os.path.join(model_path, cls._FILENAME), "rb")
            md5_file: ContextManager[TextIO] = open(os.path.join(model_path, f"{cls._FILENAME}.md5"))
        else:
            pkl_file = files(__name__).joinpath(cls._FILENAME).open("rb")
            md5_file = files(__name__).joinpath(f"{cls._FILENAME}.md5").open()
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
            pkl_file.seek(0)  # type: ignore
            return pickle.load(bin)  # type: ignore

    def __init__(
        self,
        feature_type: str = "protein",
        algorithm: str = "lbfgs",
        window_size: int = 5,
        window_step: int = 1,
        **kwargs: Any,
    ) -> None:
        """Create a new `ClusterCRF` instance.

        Arguments:
            feature_type (`str`): Defines how features should be extracted.
                Should be either *domain* or *protein*.
            algorithm (`str`): The optimization algorithm for the model. See
                https://sklearn-crfsuite.readthedocs.io/en/latest/api.html
                for available values.
            window_size (`int`): The size of the sliding window to use
                when training and predicting probabilities on sequences
                of genes.
            window_step (`int`): The step between consecutive sliding
                windows to use when training and predicting probabilities
                on sequences of genes.

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
        self.model = None
        self._options = {"algorithm": algorithm, **kwargs}

    def predict_probabilities(
        self,
        genes: Iterable[Gene],
        *,
        pad: bool = True,
        progress: Optional[Callable[[int, int], None]] = None,
    ) -> List[Gene]:
        """Predict how likely each given gene is part of a gene cluster.

        Arguments:
            genes (iterable of `~gecco.model.Gene`): The genes to compute
                probabilities for.

        Keyword Arguments:
            batch_size (`int`): The number of samples to load per batch.
                *Ignored, always 1 with the CRF.*
            pad (`bool`): Whether to pad sequences too small for a single
                window. Setting this to `False` will skip probability
                prediction entirely for sequences smaller than the window
                size.
            progress (callable): A callable that accepts two `int`, the
                current batch index and the total number of batches.

        Returns:
            `list` of `~gecco.model.Gene`: A list of new `Gene` objects with
            their probability set.

        Raises:
            `~sklearn.exceptions.NotFittedError`: When calling this method
                on an object that has not been fitted yet.

        """
        # silence progress if no callback given, ignored
        _progress = progress or (lambda x,y: None)
        window_index = 0

        # check that the model was trained
        if self.model is None:
            raise NotFittedError("This ClusterCRF instance is not fitted yet.")

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
        genes = sorted(genes, key=operator.attrgetter("source.id", "start"))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start"))

        # group genes by contigs
        contigs = {}
        for contig_id, group in itertools.groupby(genes, key=operator.attrgetter("source.id")):
            contigs[contig_id] = list(group)

        # extract features from all contigs
        contig_features = {}
        deltas = {}
        for contig_id, contig in contigs.items():
            # extract features
            feats: List[Dict[str, bool]] = extract_features(contig)
            deltas[contig_id] = 0
            # ignore sequences too small with a warning
            if len(feats) < self.window_size:
                if pad:
                    unit = self.feature_type if self.window_size - len(feats) == 1 else f"{self.feature_type}s"
                    warnings.warn(
                        f"Contig {contig[0].source.id!r} does not contain enough"
                        f" {self.feature_type}s ({len(contig)}) for sliding window"
                        f" of size {self.window_size}, padding with"
                        f" {self.window_size - len(feats)} {unit}"
                    )
                    # insert empty features as padding on both ends
                    deltas[contig_id] = delta = self.window_size - len(feats)
                    feats = [{} for _ in range(delta // 2)] + feats + [{} for _ in range((delta+1)//2)]
                else:
                    warnings.warn(
                        f"Contig {contig[0].source.id!r} does not contain enough"
                        f" {self.feature_type}s ({len(contig)}) for sliding window"
                        f" of size {self.window_size}"
                    )
                    continue
            # store features for the current contig
            contig_features[contig_id] = feats

        # compute total number of windows to process
        total = sum(len(feats) - self.window_size + 1 for feats in contig_features.values())
        _progress(window_index, total)

        # predict probabilities
        predicted = []
        for contig_id, contig in contigs.items():
            # get features if the sequence was not skipped
            if contig_id not in contig_features:
                predicted.extend(contig)
                continue
            feats = contig_features[contig_id]
            # predict marginals over a sliding window, storing maximum probabilities
            probabilities = numpy.zeros(max(len(contig), self.window_size))
            for win in sliding_window(len(feats), self.window_size, self.window_step):
                marginals = [p['1'] for p in self.model.predict_marginals_single(feats[win])]
                numpy.maximum(probabilities[win], marginals, out=probabilities[win])
                window_index += 1
                _progress(window_index, total)
            # label genes with maximal probabilities
            predicted.extend(annotate_probabilities(contig, probabilities[deltas[contig_id]//2:][:len(contig)]))

        # label domains with their biosynthetic weight according to the CRF state weights
        predicted = [
            gene.with_protein(gene.protein.with_domains(
                domain.with_cluster_weight(
                    self.model.state_features_.get((domain.name, '1'))
                )
                for domain in gene.protein.domains
            ))
            for gene in predicted
        ]

        # return the genes that were passed as input but now having 
        # gene cluster probabilities
        return predicted

    def fit(
        self,
        genes: Iterable[Gene],
        *,
        select: Optional[float] = None,
        shuffle: bool = True,
        cpus: Optional[int] = None,
        correction_method: Optional[str] = None,
    ) -> None:
        """Fit the CRF model to the given training data.

        Arguments:
            genes (iterable of `~gecco.model.Gene`): The genes to extract
                domains from for training the CRF.
            select (`float`, *optional*): The fraction of features to
                select based on Fisher-tested significance. Leave as `None`
                to skip feature selection.
            shuffle (`bool`): Whether or not to shuffle the contigs
                after having grouped the genes together.
            correction_method (`str`, *optional*): The correction method
                to use for correcting p-values used for feature selection.
                Ignored if ``select`` is `False`.

        """
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
            sorted_sig = sorted(sig, key=sig.get)[:int(select*len(sig))]  # type: ignore
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
        self.model = model = sklearn_crfsuite.CRF(**self._options)
        model.fit(training_features, training_labels)

    def save(self, model_path: "os.PathLike[str]") -> None:
        """Save the `ClusterCRF` to an on-disk location.

        Models serialized at a given location can be later loaded from that
        same location using the `ClusterCRF.trained` class method.

        Arguments:
            model_path (`str`): The path to the directory where to write
                the model files.

        """
        # pickle the CRF model
        model_out = os.path.join(model_path, self._FILENAME)
        with open(model_out, "wb") as out:
            pickle.dump(self, out, protocol=4)
        # compute a checksum for the pickled model
        hasher = hashlib.md5()
        with open(model_out, "rb") as out:
            for chunk in iter(lambda: out.read(io.DEFAULT_BUFFER_SIZE), b""):
                hasher.update(chunk)
        # write the checksum next to the pickled file
        with open(f"{model_out}.md5", "w") as out_hash:
            out_hash.write(hasher.hexdigest())
