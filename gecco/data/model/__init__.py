import functools
import hashlib
import io
import pickle
import pkg_resources
import typing

if typing.TYPE_CHECKING:
    from ...crf import ClusterCRF


def load() -> 'ClusterCRF':
    """Unpickle the CRF model from `gecco.data` after checking its hash.
    """
    with pkg_resources.resource_stream(__name__, "crf.model.md5") as sig:
        signature = sig.read().decode("ascii").strip()
    hasher = hashlib.md5()
    with pkg_resources.resource_stream(__name__, "crf.model") as bin:
        read = functools.partial(bin.read, io.DEFAULT_BUFFER_SIZE)
        for chunk in iter(read, b""):
            hasher.update(typing.cast(bytes, chunk))
    if hasher.hexdigest().upper() != signature.upper():
        raise RuntimeError("hashes of model data does not match signature")
    with pkg_resources.resource_stream(__name__, "crf.model") as bin:
        return pickle.load(bin)  # type: ignore
