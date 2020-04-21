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
            self.logger._log(verboselogs.NOTICE, self.Message(msg, args), **kw)

    def spam(self, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(verboselogs.SPAM):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(verboselogs.SPAM, self.Message(msg, args), **kw)

    def verbose(self, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(verboselogs.VERBOSE):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(verboselogs.VERBOSE, self.Message(msg, args), **kw)

    def success(self, msg: str, *args: object, **kwargs: Any) -> None:
        if self.isEnabledFor(verboselogs.SUCCESS):
            msg, kw = self.process(msg, kwargs)
            self.logger._log(verboselogs.SUCCESS, self.Message(msg, args), **kw)


def wrap_warnings(logger: logging.Logger) -> Callable[["_F"], "_F"]:
    """Have the function patch `warnings.showwarning` with the given logger.

    Arguments:
        logger (~logging.Logger): the logger to wrap warnings with when
            the decorated function is called.

    Returns:
        `function`: a decorator function that will wrap a callable and
        redirect any warning raised by that callable to the given logger.

    Example:
        >>> logger = logging.getLogger()
        >>> @wrap_warnings(logger)
        ... def divide_by_zero(x):
        ...     return numpy.array(x) / 0

    """

    class _WarningsWrapper(object):
        def __init__(self, logger: logging.Logger, func: Callable[..., "_T"]):
            self.logger = logger
            self.func = func
            functools.update_wrapper(self, func)

        def showwarning(
            self,
            message: str,
            category: Type[Warning],
            filename: str,
            lineno: int,
            file: Optional[TextIO] = None,
            line: Optional[str] = None,
        ) -> None:
            for line in filter(str.strip, str(message).splitlines()):
                self.logger.warning(line.strip())

        def __call__(self, *args: Any, **kwargs: Any) -> "_T":
            old_showwarning = warnings.showwarning
            warnings.showwarning = self.showwarning
            try:
                return self.func(*args, **kwargs)
            finally:
                warnings.showwarning = old_showwarning

        def __getattr__(self, name: Any) -> Any:
            return getattr(self.func, name)

    def decorator(func: Callable[..., "_T"]) -> Callable[..., "_T"]:
        return _WarningsWrapper(logger, func)

    return decorator  # type: ignore


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
