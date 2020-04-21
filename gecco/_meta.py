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
    """

    def __init__(self, f: Callable[["_S"], "_T"]) -> None:
        self.f = f

    def __get__(self, obj: object, owner: "_S") -> "_T":  # type: ignore
        return self.f(owner)
