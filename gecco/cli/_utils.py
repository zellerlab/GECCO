"""A module containing some utilities only relevant to the CLI.
"""

import contextlib
import functools
import logging
import typing
import warnings
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Type, TextIO

import numpy
import verboselogs

from .._meta import classproperty

if typing.TYPE_CHECKING:
    _S = typing.TypeVar("_S")
    _T = typing.TypeVar("_T")
    _F = typing.TypeVar("_F", bound=Callable[..., "_T"])
    ShowWarning = typing.Callable[
        [str, Type[Warning], str, int, Optional[TextIO], Optional[str]],
        None
    ]


class BraceAdapter(logging.LoggerAdapter, verboselogs.VerboseLogger):
    """A logging adapter for `VerboseLogger` to use new-style formatting.
    """

    class Message(object):
        def __init__(self, fmt: object, args: Iterable[object]):
            self.fmt = str(fmt)
            self.args = args

        def __str__(self) -> str:
            return self.fmt.format(*self.args)

    def __init__(
        self, logger: logging.Logger, extra: Optional[Dict[str, object]] = None
    ) -> None:
        super(BraceAdapter, self).__init__(logger, extra or {})
        self.handlers = self.logger.handlers

    @property
    def level(self) -> int:
        return self.logger.level

    def log(self, level: int, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(level):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(level, self.Message(msg, args), (), **kw)

    def notice(self, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(verboselogs.NOTICE):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(verboselogs.NOTICE, self.Message(msg, args), (), **kw)

    def spam(self, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(verboselogs.SPAM):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(verboselogs.SPAM, self.Message(msg, args), (), **kw)

    def verbose(self, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(verboselogs.VERBOSE):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(verboselogs.VERBOSE, self.Message(msg, args), (), **kw)

    def success(self, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(verboselogs.SUCCESS):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(verboselogs.SUCCESS, self.Message(msg, args), (), **kw)


@contextlib.contextmanager
def patch_showwarnings(new_showwarning: "ShowWarning") -> Iterator[None]:
    """Make a context patching `warnings.showwarning` with the given function.
    """
    old_showwarning = warnings.showwarning
    try:
        warnings.showwarning = new_showwarning
        yield
    finally:
        warnings.showwarning = old_showwarning


@contextlib.contextmanager
def numpy_error_context(
    *,
    all: Optional[str] = None,
    divide: Optional[str] = None,
    over: Optional[str] = None,
    under: Optional[str] = None,
    invalid: Optional[str] = None
) -> Iterator[None]:
    """A context manager to modify the `numpy` error behaviour locally.

    Example:
        >>> with numpy_error_context(divide="ignore"):
        ...     numpy.log10(0)
        -inf
        >>> with numpy_error_context(divide="raise"):
        ...     numpy.log10(0)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        FloatingPointError: divide by zero encountered in log10

    """
    try:
        old_settings = numpy.seterr(
            all=all, divide=divide, over=over, under=under, invalid=invalid
        )
        yield
    finally:
        numpy.seterr(**old_settings)


def guess_sequences_format(path: str) -> Optional[str]:
    """Guess the format of a sequence file located in ``path``.

    Supports the following formats:
        * GenBank
        * FASTA

    """
    head = ""
    with open(path, "r") as file:
        for head in iter(lambda: file.read(256), ""):
            head = head.strip()
            if head:
                break
    if head.startswith(">"):
        return "fasta"
    elif head.startswith("LOCUS"):
        return "genbank"
    else:
        return None


def in_context(func):

    @functools.wraps(func)
    def newfunc(*args, **kwargs):
        with contextlib.ExitStack() as ctx:
            func(*args, ctx, **kwargs)

    return newfunc
