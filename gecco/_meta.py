"""Share metaprogramming helpers for GECCO.
"""

import abc
import typing
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
