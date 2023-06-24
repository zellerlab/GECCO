"""A module containing some utilities only relevant to the CLI.
"""

import contextlib
import functools
import io
import logging
import typing
import warnings
from types import ModuleType, TracebackType
from typing import Any, BinaryIO, Callable, Dict, Iterable, Iterator, List, Optional, Type, TextIO, Union, Tuple

from .._meta import classproperty, zopen

if typing.TYPE_CHECKING:
    from rich.progress import Progress, TaskID

    ShowWarning = typing.Callable[
        [Union[Warning, str], Type[Warning], str, int, Optional[TextIO], Optional[str]],
        None
    ]
    _F = typing.TypeVar("_F", bound=Callable[..., Any])


class ProgressReader(io.RawIOBase):
    """A reader that updates a progress bar while it's being read from.
    """

    @staticmethod
    def scale_size(length: float) -> Tuple[float, int, str]:
        for scale, unit in enumerate(["B", "kiB", "MiB", "GiB", "TiB"]):
            if length > 1024:
                length /= 1024
            else:
                break
        return length, scale, unit

    def __init__(
        self,
        handle: BinaryIO,
        progress: "Progress",
        task: "TaskID",
        scale: int = 0
    ):
        self.handle = handle
        self.progress = progress
        self.task = task
        self.scale = scale

    def __enter__(self) -> "ProgressReader":
        self.handle.__enter__()
        return self

    def __exit__(
        self,
        exc_ty: Optional[Type[BaseException]],
        exc_val: Optional[BaseException],
        tb: Optional[TracebackType]
    ) -> None:
        self.handle.__exit__(exc_ty, exc_val, tb)

    def _update(self, length: int) -> None:
        self.progress.update(self.task, advance=length / (1024 ** self.scale))

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return False

    def writable(self) -> bool:
        return False

    def readline(self, size: Optional[int] = -1) -> bytes:
        line = self.handle.readline(-1 if size is None else size)
        self._update(len(line))
        return line

    def readlines(self, hint: Optional[int] = -1) -> List[bytes]:
        lines = self.handle.readlines(-1 if hint is None else hint)
        self._update(sum(map(len, lines)))
        return lines

    def read(self, size: Optional[int] = -1) -> bytes:
        block = self.handle.read(-1 if size is None else size)
        self._update(len(block))
        return block

    def close(self) -> None:
        self.handle.close()


@contextlib.contextmanager
def patch_showwarnings(new_showwarning: "ShowWarning") -> Iterator[None]:
    """Make a context patching `warnings.showwarning` with the given function.
    """
    old_showwarning: ShowWarning = warnings.showwarning
    try:
        warnings.showwarning = new_showwarning  # type: ignore
        yield
    finally:
        warnings.showwarning = old_showwarning  # type: ignore


@contextlib.contextmanager
def numpy_error_context(
    numpy: ModuleType,
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
        * EMBL

    """
    head = ""
    with zopen(path) as file:
        for head in iter(lambda: file.read(256), b""):
            head = head.strip()
            if head:
                break
    if head.startswith(b">"):
        return "fasta"
    elif head.startswith(b"LOCUS"):
        return "genbank"
    elif head.startswith(b"ID "):
        return "embl"
    else:
        return None


def in_context(func: "_F") -> "_F":

    @functools.wraps(func)
    def newfunc(*args, **kwargs):  # type: ignore
        with contextlib.ExitStack() as ctx:
            return func(*args, ctx, **kwargs)

    return newfunc  # type: ignore
