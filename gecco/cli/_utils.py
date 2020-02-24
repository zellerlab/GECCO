import contextlib
import typing

import numpy


@contextlib.contextmanager
def numpy_error_context(**kwargs: typing.Dict[str, str]):
    """A context manager to modify the `numpy` error behaviour locally.

    Example:
        >>> with numpy_error_context(divide="ignore"):
        ...     numpy.log10(0)
        -inf
        >>> with numpy_error_context(divide="error"):
        ...     numpy.log10(0)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        FloatingPointError: divide by zero encountered in log10

    """
    try:
        old_settings = numpy.seterr(**kwargs)
        yield
    finally:
        numpy.seterr(**old_settings)
