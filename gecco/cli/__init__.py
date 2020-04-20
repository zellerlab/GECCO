# coding: utf-8

import sys
import typing
from typing import Optional, List, TextIO

import demandimport

if typing.TYPE_CHECKING:
    import logging

def main(
    argv: Optional[List[str]] = None,
    stream: Optional[TextIO] = None,
    logger: Optional["logging.Logger"] = None,
) -> int:
    # enable demandimport only when importing the command and parsing the
    # arguments, but disable it for actual execution of the app
    with demandimport.enabled():
        demandimport.ignore('msvcrt')
        demandimport.ignore('_compat_pickle')
        from .commands._main import Main
        _main = Main(argv, stream, logger)
    return _main()
