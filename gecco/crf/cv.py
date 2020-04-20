import math
import multiprocessing
from typing import Iterable, Iterator, List, Tuple

import numpy
import sklearn.model_selection


class LotoSplit(object):

    def __init__(self, type_array: Iterable[List[str]]) -> None:
        self.type_array = list(type_array)
        self.type_set = {ty for types in self.type_array for ty in types}

    def split(self) -> Iterator[Tuple["numpy.ndarray", "numpy.ndarray", str]]:
        for typ in self.type_set:
            train_idx = numpy.array([idx for idx, s in enumerate(self.type_array)
                if typ not in s])
            test_idx = numpy.array([idx for idx, s in enumerate(self.type_array)
                if typ in s])
            yield train_idx, test_idx, typ


class StratifiedSplit(object):

    def __init__(self, type_array: Iterable[List[str]], n_splits: int=10) -> None:
        self.type_array = [t[0] if len(t) == 1 else "Mixed" for t in type_array]
        type_set = list(set(self.type_array))
        self.type_enc = [type_set.index(t) for t in self.type_array]
        self.skf = sklearn.model_selection.StratifiedKFold(n_splits=n_splits)

    def split(self) -> Iterator[Tuple["numpy.ndarray", "numpy.ndarray"]]:
        return self.skf.split(self.type_array, self.type_enc)  # type: ignore


# FUNC
def n_folds(n_samples, n=10):
    """Generates test fold for sklearn's PredefinedSplit."""
    fold_size = math.ceil(n_samples / n)
    indices = numpy.array(range(0, math.ceil(n_samples / fold_size)))
    folds = numpy.repeat(indices, fold_size)[:n_samples]
    return folds


def n_folds_partial(n_train, n_val, n=10):
    """
    Generates test fold for sklearn's PredefinedSplit so that
    n_train are only used for training and
    n_val are used to do cross_validation.
    Creates <n_folds> folds from the cv_set.
    """
    train_folds = numpy.array([-1] * n_train)
    fold_size = math.ceil(n_val / n)
    cv_indices = numpy.array(range(0, math.ceil(n_val / fold_size)))
    cv_folds = numpy.repeat(cv_indices, fold_size)[:n_val]
    return numpy.append(cv_folds, train_folds)


def single_fold(train_set, cv_set, fold, n=10):
    """Creates a single fold from n_folds."""
    folds = n_folds(train_set, cv_set, n)
    folds[folds != fold] = -1
    folds[folds == fold] = 1
    return folds


def fold_mask(train_set, cv_set, fold, n=10):
    """Creates a fold mask to select train and test set from array."""
    folds = n_folds(train_set, cv_set, n)
    folds[folds != fold] = True
    folds[folds == fold] = False
    return folds


def leave_n_out_folds(train_set, cv_set, n=1):
    """
    Generates test fold for sklearn's PredefinedSplit so that
    train_set is only used for training and
    cv_set is used to do cross_validation.
    n determines how many samples are left out at each iteration.
    """
    train_folds = numpy.array([-1] * len(train_set))
    cv_indices = numpy.array(range(0, math.ceil(len(cv_set) / n)))
    cv_folds = numpy.repeat(cv_indices, n)[:len(cv_set)]
    return numpy.append(train_folds, cv_folds)
