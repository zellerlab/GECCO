import io
import os
import pkg_resources
import typing
from typing import BinaryIO, IO, Union

from . import hmms, model, knn


def realpath(local_path: str) -> str:
    """Get the system path to a data file in the `gecco.data` module.

    Raises:
        FileNotFoundError: when the data file cannot be located.
    """
    path = pkg_resources.resource_filename(__name__, local_path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return path


def open(local_path: str, mode: str = "r") -> Union[IO[str], IO[bytes]]:
    """Get a file handle to a data file in the `gecco.data` module.
    """
    stream = pkg_resources.resource_stream(__name__, local_path)
    if mode == "r":
        stream = io.TextIOWrapper(stream)  # type: ignore
    elif mode not in {"rb", "br"}:
        raise ValueError(f"invalid mode: {mode!r}")
    return stream
