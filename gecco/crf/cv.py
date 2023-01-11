"""Cross-validation utilities for `~gecco.crf.ClusterCRF`.
"""

import math
import multiprocessing
import typing
from typing import Any, Iterable, Iterator, List, Set, Tuple

import numpy
import sklearn.model_selection

if typing.TYPE_CHECKING:
    from numpy.typing import NDArray


class LeaveOneGroupOut(sklearn.model_selection.LeaveOneGroupOut):
    """A `~sklearn.model_selection.LeaveOneGroupOut` supporting multiple labels.

    If a sample has multiple class labels, it will be excluded from both
    training and testing data when one of its labels corresponds to the
    fold::

        >>> loto = LeaveOneGroupOut()
        >>> groups = [["a"], ["b"], ["c"], ["a", "b"]]
        >>> for i, (trn, tst) in enumerate(loto.split(range(4), groups=groups)):
        ...     print("-"*20)
        ...     print(" FOLD", i+1)
        ...     print("TRAIN", f"{str(trn):<7}", [groups[i] for i in trn])
        ...     print(" TEST", f"{str(tst):<7}", [groups[i] for i in tst])
        ...
        --------------------
         FOLD 1
        TRAIN [1 2]   [['b'], ['c']]
         TEST [0]     [['a']]
        --------------------
         FOLD 2
        TRAIN [0 2]   [['a'], ['c']]
         TEST [1]     [['b']]
        --------------------
         FOLD 3
        TRAIN [0 1 3] [['a'], ['b'], ['a', 'b']]
         TEST [2]     [['c']]

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

    def split(self, X: Any, y: Any = None, groups: Any = None) -> Iterator[Tuple["NDArray[numpy.int_]", "NDArray[numpy.int_]"]]:  # noqa: D102
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None")
        # collect groups
        group_sets: List[Set[object]] = list(map(set, groups))
        unique_groups = {label for labels in group_sets for label in labels}
        #
        indices = numpy.arange(len(X))
        for ty in sorted(unique_groups):  # type: ignore
            test_mask = numpy.array([list(group) == [ty] for group in groups])
            train_mask = numpy.array([ty not in group for group in groups])
            yield indices[train_mask], indices[test_mask]
