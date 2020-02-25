import hashlib
import io
import os
import pickle
import pkg_resources
import typing


def realpath(local_path: str) -> str:
    """Get the system path to a data file in the `gecco.data` module.

    Raises:
        FileNotFoundError: when the data file cannot be located.
    """
    path = pkg_resources.resource_filename(__name__, local_path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return path


def open(local_path: str, mode: str = "r") -> typing.TextIO:
    """Get a file handle to a data file in the `gecco.data` module.
    """
    stream = pkg_resources.resource_stream(__name__, local_path)
    if mode == "r":
        stream = io.TextIOWrapper(stream)
    elif mode not in {"rb", "br"}:
        raise ValueError(f"invalid mode: {mode!r}")
    return stream


def load(local_path: str) -> object:
    """Unpickle an object from the `gecco.data` after checking its hash.
    """
    hasher = hashlib.md5()
    with open(f"{local_path}.md5") as sig:
        signature = sig.read().strip()
    with open(local_path, "rb") as bin:
        read = lambda: bin.read(io.DEFAULT_BUFFER_SIZE)
        for chunk in iter(read, b''):
            hasher.update(chunk)
    if hasher.hexdigest().upper() != signature.upper():
        raise RuntimeError("hashes of data does not match signature")
    with open(local_path, "rb") as bin:
        return pickle.load(bin)
