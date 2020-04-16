# coding: utf-8

import sys
import typing
from typing import Optional, List, TextIO

import verboselogs
import lazy_import

# Add here all the external modules that could be slow to load
lazy_import.lazy_module("Bio.SeqIO")
lazy_import.lazy_module("numpy")
lazy_import.lazy_module("pandas")
lazy_import.lazy_module("sklearn_crfsuite")
lazy_import.lazy_module("sklearn.model_selection")
lazy_import.lazy_module("sklearn.neighbors")
lazy_import.lazy_module("scipy.spatial.distance")
lazy_import.lazy_module("tqdm")

from .commands._main import Main as _Main

if typing.TYPE_CHECKING:
    import logging


def main(
    argv: Optional[List[str]] = None,
    stream: Optional[TextIO] = None,
    logger: Optional["logging.Logger"] = None,
) -> int:
    return _Main(argv, stream, logger)()

if __name__ == "__main__":
    sys.exit(main())
