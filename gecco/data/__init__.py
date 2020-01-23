import io
import pkg_resources


def realpath(local_path):
    return pkg_resources.resource_filename(__name__, local_path)


def open(local_path, mode="r"):
    stream = pkg_resources.resource_stream(__name__, local_path)
    if mode == "r":
        stream = io.TextIOWrapper(stream)
    elif mode not in {"rb", "br"}:
        raise ValueError(f"invalid mode: {mode!r}")
    return stream
