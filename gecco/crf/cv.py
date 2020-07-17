"""Cross-validation utilities for `~gecco.crf.ClusterCRF`.
"""

import math
import multiprocessing
from typing import Iterable, Iterator, List, Set, Tuple

import numpy
import sklearn.model_selection


class LeaveOneGroupOut(sklearn.model_selection.LeaveOneGroupOut):
    """A `~sklearn.model_selection.LeaveOneGroupOut` supporting multiple labels.

    If a sample has multiple class labels, it will be left out of the training
    data each time the fold corresponds to one of its class labels::

        >>> loto = LeaveOneGroupOut()
        >>> groups = numpy.array([ ["a"], ["b"], ["c"], ["a", "b"] ])
        >>> for i, (trn, tst) in enumerate(loto.split(range(4), groups=groups)):
        ...     print("-"*20)
        ...     print(" FOLD", i+1)
        ...     print("TRAIN", trn, list(groups[trn]))
        ...     print(" TEST", tst, list(groups[tst]))
        ...
        --------------------
         FOLD 1
        TRAIN [1 2] [['b'], ['c']]
         TEST [0 3] [['a'], ['a', 'b']]
        --------------------
         FOLD 2
        TRAIN [0 2] [['a'], ['c']]
         TEST [1 3] [['b'], ['a', 'b']]
        --------------------
         FOLD 3
        TRAIN [0 1 3] [['a'], ['b'], ['a', 'b']]
         TEST [2] [['c']]

    """

    def get_n_splits(
        self, X: object = None, y: object = None, groups: Iterable[Iterable[str]] = None
    ) -> int:
        """Return the number of splitting iterations in the cross-validator.

        Arguments:
            X (object): Always ignored, exists for compatibility.
            Y (object): Always ignored,
            groups (iterable of iterable of object): Group labels for the
                samples used while splitting the dataset into train/test set.
                It must contain `n_samples` elements each having one or more
                labels. This ``groups`` parameter must always be specified to
                calculate the number of splits, though the other parameters
                can be omitted.

        Returns:
            int: The number of splitting iterations in the cross-validator.
            *This is equal to the number of unique classes within the groups*.

        Raises:
            ValueError: When ``groups`` is not given.

        Example:
            >>> loto = LeaveOneGroupOut()
            >>> groups = [["Polyketide"], ["NRP"], ["RiPP"]]
            >>> loto.get_n_splits(groups=groups)
            3
            >>> groups = [["Terpene"], ["NRP"], ["RiPP"], ["Terpene", "NRP"]]
            >>> loto.get_n_splits(groups=groups)
            3

        """
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None")
        labels = {label for labels in groups for label in labels}
        return len(labels)

    def _iter_test_masks(
        self, X: object, y: object, groups: Iterable[Iterable[object]]
    ) -> Iterator["numpy.ndarray"]:
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        # We collect the groups to avoid side-effects during iteration
        group_sets: List[Set[object]] = list(map(set, groups))  # type: ignore
        unique_groups = {label for labels in group_sets for label in labels}
        if len(unique_groups) <= 1:
            raise ValueError(
                "The groups parameter contains fewer than 2 unique groups "
                f"({unique_groups}). LeaveOneGroupOut expects at least 2."
            )
        for i in sorted(unique_groups):
            yield numpy.array([i in group for group in groups])
