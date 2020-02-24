import io
import os
import pkg_resources


def realpath(local_path):
    """Get the system path to a data file in the `gecco.data` module.

    Raises:
        FileNotFoundError: when the data file cannot be located.
    """
    path = pkg_resources.resource_filename(__name__, local_path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return path


def open(local_path, mode="r"):
    """Get a file handle to a data file in the `gecco.data` module.
    """
    stream = pkg_resources.resource_stream(__name__, local_path)
    if mode == "r":
        stream = io.TextIOWrapper(stream)
    elif mode not in {"rb", "br"}:
        raise ValueError(f"invalid mode: {mode!r}")
    return stream
