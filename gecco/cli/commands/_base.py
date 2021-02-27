# coding: utf-8
import abc
import datetime
import logging
import os
import socket
import sys
import textwrap
import typing
from typing import Any, ClassVar, Callable, Optional, List, Mapping, Dict, TextIO, Type

import docopt
import rich.console
import rich.logging

from ... import __version__, __name__ as __progname__
from ._error import CommandExit, InvalidArgument

_T = typing.TypeVar("_T")

class Command(metaclass=abc.ABCMeta):
    """An abstract base class for ``gecco`` subcommands.
    """

    # -- Abstract methods ----------------------------------------------------

    doc: ClassVar[str] = NotImplemented
    summary: ClassVar[str] = NotImplemented

    @abc.abstractmethod
    def __call__(self) -> int:
        return NotImplemented  # type: ignore

    @classmethod
    @abc.abstractmethod
    def doc(cls, fast: bool = False) -> str:
        return NotImplemented  # type: ignore

    # -- Concrete methods ----------------------------------------------------

    _version: str = "{} {}".format(__progname__, __version__)
    _options_first: bool = False

    def __init__(
        self,
        argv: Optional[List[str]] = None,
        stream: Optional[TextIO] = None,
        options: Optional[Mapping[str, Any]] = None,
        config: Optional[Dict[Any, Any]] = None,
    ) -> None:
        self._stream: Optional[TextIO] = stream
        self.argv = argv
        self.stream: TextIO = stream or sys.stderr
        self.options = options or dict()
        self.pool = None
        self.config = config
        self.console = rich.console.Console(file=self.stream)

        self._hostname = socket.gethostname()
        self._pid = os.getpid()

        # Parse command line arguments
        try:
            self.args = docopt.docopt(
                textwrap.dedent(self.doc(fast=True)).lstrip(),
                help=False,
                argv=argv,
                version=self._version,
                options_first=self._options_first,
            )
            self.verbose = self.args.get("--verbose", 0)
            self.quiet = self.args.get("--quiet", 0)
        except docopt.DocoptExit as de:
            self.args = de
            self.level = 0
            self.quiet = 0

    def _check(self) -> None:
        # Assert CLI arguments were parsed Successfully
        if isinstance(self.args, docopt.DocoptExit):
            self.console.print(self.args)#, file=self.stream)
            raise CommandExit(1)
        return None

    def _check_flag(
        self,
        name: str,
        convert: Optional[Callable[[str], _T]] = None,
        check: Optional[Callable[[_T], bool]] = None,
        message: Optional[str] = None,
        hint: Optional[str] = None,
    ) -> _T:
        _convert = (lambda x: x) if convert is None else convert
        _check = (lambda x: True) if check is None else check
        try:
            value = _convert(self.args[name])
            if not _check(value):
                raise ValueError(self.args[name])
        except Exception as err:
            if hint is None:
                self.error(f"Invalid value for argument [purple]{name}[/]:", repr(self.args[name]))
            else:
                self.error(f"Invalid value for argument [purple]{name}[/]:", repr(self.args[name]), f"(expected {hint})")
            raise InvalidArgument(self.args[name]) from err
        else:
            return value

    # -- Logging methods -----------------------------------------------------

    def error(self, message, *args, level=0):
        if self.quiet <= 2 and level <= self.verbose:
            self.console.print(
                *self._logprefix(),
                "[bold red]FAIL[/]",
                message,
                *args,
            )

    def info(self, verb, *args, level=1):
        if self.quiet == 0 and level <= self.verbose:
            self.console.print(
                *self._logprefix(),
                f"[bold blue]INFO[/]",
                verb,
                *args,
            )

    def success(self, verb, *args, level=1):
        if self.quiet == 0 and level <= self.verbose:
            self.console.print(
                *self._logprefix(),
                f"[bold green]  OK[/]",
                verb,
                *args,
            )

    def warn(self, verb, *args, level=0):
        if self.quiet <= 1 and level <= self.verbose:
            self.console.print(
                *self._logprefix(),
                "[bold yellow]WARN[/]",
                verb,
                *args
            )

    def _logprefix(self):
        return [
            f"[dim cyan]{datetime.datetime.now().strftime('%Y-%m-%d %H:%m:%S')}[/]",
            f"[dim purple]{self._hostname}[/]",
            f"[dim]{__progname__}[[default dim]{self._pid}[/]][/]",
        ]

    def _showwarnings(
        self,
        message: str,
        category: Type[Warning],
        filename: str,
        lineno: int,
        file: Optional[TextIO] = None,
        line: Optional[str] = None,
    ) -> None:
        for line in filter(str.strip, str(message).splitlines()):
            self.warn(line.strip())
