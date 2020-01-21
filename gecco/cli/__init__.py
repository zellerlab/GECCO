# coding: utf-8

import logging
import sys
import typing

import verboselogs

from .commands._main import Main as _Main


def main(
    argv: typing.Optional[typing.List[str]] = None,
    stream: typing.Optional[typing.TextIO] = None,
    logger: typing.Optional[logging.Logger] = None,
) -> int:
    return _Main(argv, stream, logger)()


if __name__ == "__main__":
    sys.exit(main())
