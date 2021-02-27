# coding: utf-8

import sys
import typing
from typing import Optional, List, TextIO

if typing.TYPE_CHECKING:
    import logging


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:
    from .commands._main import Main

    return Main(argv, stream).execute()
