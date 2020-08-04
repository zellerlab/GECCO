"""Shared metaprogramming helpers for GECCO.
"""

import abc
import functools
import operator
import typing
from multiprocessing.pool import Pool
from typing import Callable

if typing.TYPE_CHECKING:
    _S = typing.TypeVar("_S")
    _T = typing.TypeVar("_T")


class classproperty(property):
    """A class property decorator.

    Example:
        Use `classproperty` to create a property which counts the number of
        times it has been accessed:

        >>> class X():
        ...     __COUNT = 0
        ...     @classproperty
        ...     def count(cls):
        ...         cls.__COUNT += 1
        ...         return cls.__COUNT
        >>> X.count
        1
        >>> X.count
        2

    """

    def __init__(self, f: Callable[["_S"], "_T"]) -> None:
        self.f = f

    def __get__(self, obj: object, owner: "_S") -> "_T":  # type: ignore
        return self.f(owner)


class OrderedPoolWrapper:
    """A `Pool` wrapper that returns results in the order they were given.
    """

    class _OrderedFunc:

        def __init__(self, inner, star=False):
            self.inner = inner
            self.star = star

        def __call__(self, args):
            i, other = args
            if self.star:
                return i, self.inner(*other)
            else:
                return i, self.inner(other)

    def __init__(self, inner: Pool) -> None:
        self.inner = inner

    def __enter__(self) -> "OrderedPool":
        self.inner.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        return self.inner.__exit__(exc_type, exc_value, exc_tb)

    def map(self, func, it):
        wrapped_it = enumerate(it)
        wrapped_func = self._OrderedFunc(func)
        results = self.inner.map(wrapped_func, wrapped_it)
        results.sort(key=operator.itemgetter(0))
        return list(map(operator.itemgetter(1), results))

    def starmap(self, func, it):
        wrapped_it = enumerate(it)
        wrapped_func = self._OrderedFunc(func, star=True)
        results = self.inner.map(wrapped_func, wrapped_it)
        results.sort(key=operator.itemgetter(0))
        return list(map(operator.itemgetter(1), results))
