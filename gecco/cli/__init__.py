# coding: utf-8

import contextlib
import sys
import typing
from typing import Optional, List, TextIO

def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:
    from .commands._main import Main

    with contextlib.ExitStack() as ctx:
        return Main(argv, stream).execute(ctx)

    return 0
