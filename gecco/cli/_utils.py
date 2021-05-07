"""A module containing some utilities only relevant to the CLI.
"""

import contextlib
import functools
import logging
import typing
import warnings
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Type, TextIO

from .._meta import classproperty

if typing.TYPE_CHECKING:
    ShowWarning = typing.Callable[
        [str, Type[Warning], str, int, Optional[TextIO], Optional[str]],
        None
    ]

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
    numpy,
    *,
    all: Optional[str] = None,
    divide: Optional[str] = None,
    over: Optional[str] = None,
    under: Optional[str] = None,
    invalid: Optional[str] = None
) -> Iterator[None]:
    """A context manager to modify the `numpy` error behaviour locally.

    Example:
        >>> import numpy
        >>> with numpy_error_context(numpy, divide="ignore"):
        ...     numpy.log10(0)
        -inf
        >>> with numpy_error_context(numpy, divide="raise"):
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
            return func(*args, ctx, **kwargs)

    return newfunc
