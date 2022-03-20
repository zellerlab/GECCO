"""Shared metaprogramming helpers for GECCO.
"""

import abc
import contextlib
import functools
import importlib
import locale
import operator
import typing
from multiprocessing.pool import Pool
from typing import Any, Callable, Iterable, Iterator, List, Tuple, Optional, Type

if typing.TYPE_CHECKING:
    from types import TracebackType

    _S = typing.TypeVar("_S")
    _T = typing.TypeVar("_T")
    _A = typing.TypeVar("_A")
    _R = typing.TypeVar("_R")
    # _F = typing.TypeVar("_F", bound=Callable[[_A], _R])

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


class requires:
    """A decorator for functions that require optional dependencies.
    """

    def __init__(self, module_name):
        self.module_name = module_name

        try:
            self.module = importlib.import_module(module_name)
        except ImportError as err:
            self.module = err

    def __call__(self, func):

        if isinstance(self.module, ImportError):

            @functools.wraps(func)
            def newfunc(*args, **kwargs):
                msg = f"calling {func.__qualname__} requires module {self.module.name}"
                raise RuntimeError(msg) from self.module
        else:

            newfunc = func
            basename = self.module_name.split(".")[-1]
            newfunc.__globals__[basename] = self.module

        return newfunc


class UniversalContainer(object):
    """A container that contains everything.
    """

    def __repr__(self):
        return f"{self.__class__.__name__}()"

    def __contains__(self, item: object) -> bool:
        return True


def sliding_window(length: int, window: int, step: int) -> Iterator[slice]:
    """Iterate over a sequence of length `length` with a sliding window.
    """
    if window <= 0:
        raise ValueError("Window size must be strictly positive")
    if step <= 0 or step > window:
        raise ValueError("Window step must be strictly positive and under `window_size`")
    for i in range(0, length + 1 - window, step):
        yield slice(i, i+window)


@contextlib.contextmanager
def patch_locale(name: str):
    """Create a context manager to locally change the locale in use.
    """
    lc = locale.setlocale(locale.LC_TIME)
    try:
        locale.setlocale(locale.LC_TIME, name)
        yield
    finally:
        locale.setlocale(locale.LC_TIME, lc)
