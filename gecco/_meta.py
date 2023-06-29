"""Shared metaprogramming helpers for GECCO.
"""

import abc
import bz2
import contextlib
import functools
import importlib
import locale
import lzma
import operator
import typing
from multiprocessing.pool import Pool
from typing import (
    Any,
    BinaryIO,
    Callable,
    Iterable,
    Iterator,
    List,
    Tuple,
    Union,
    Optional,
    Type
)

try:
    from isal import igzip as gzip
except ImportError:
    import gzip

try:
    import lz4.frame
except ImportError as err:
    lz4 = err


if typing.TYPE_CHECKING:
    from types import TracebackType, ModuleType

    _S = typing.TypeVar("_S")
    _T = typing.TypeVar("_T")
    _A = typing.TypeVar("_A")
    _R = typing.TypeVar("_R")
    _F = typing.TypeVar("_F", bound=Callable[..., Any])

_BZ2_MAGIC = b"BZh"
_GZIP_MAGIC = b"\x1f\x8b"
_XZ_MAGIC = b"\xfd7zXZ"
_LZ4_MAGIC = b"\x04\x22\x4d\x18"


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

    module: Union["ModuleType", BaseException]

    def __init__(self, module_name: str) -> None:
        self.module_name = module_name

        try:
            self.module = importlib.import_module(module_name)
        except ImportError as err:
            self.module = err

    def __call__(self, func: "_F") -> "_F":

        if isinstance(self.module, ImportError):

            @functools.wraps(func)
            def newfunc(*args, **kwargs):  # type: ignore
                msg = f"calling {func.__qualname__} requires module {self.module.name}"
                raise RuntimeError(msg) from self.module
        else:

            newfunc = func
            basename = self.module_name.split(".")[-1]
            newfunc.__globals__[basename] = self.module

        return newfunc # type: ignore


class UniversalContainer(object):
    """A container that contains everything.
    """

    def __repr__(self) -> str:
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
def patch_locale(name: str) -> Iterator[None]:
    """Create a context manager to locally change the locale in use.
    """
    lc = locale.setlocale(locale.LC_TIME)
    try:
        locale.setlocale(locale.LC_TIME, name)
        yield
    finally:
        locale.setlocale(locale.LC_TIME, lc)


@contextlib.contextmanager
def zopen(path: str) -> Iterator[BinaryIO]:
    with contextlib.ExitStack() as ctx:
        file = ctx.enter_context(open(path, "rb"))
        peek = file.peek()
        if peek.startswith(_GZIP_MAGIC):
            file = ctx.enter_context(gzip.open(file, mode="rb"))
        elif peek.startswith(_BZ2_MAGIC):
            file = ctx.enter_context(bz2.open(file, mode="rb"))
        elif peek.startswith(_XZ_MAGIC):
            file = ctx.enter_context(lzma.open(file, mode="rb"))
        elif peek.startswith(_LZ4_MAGIC):
            if isinstance(lz4, ImportError):
                raise RuntimeError("File compression is LZ4 but python-lz4 is not installed") from lz4
            file = ctx.enter_context(lz4.frame.open(file))
        yield file