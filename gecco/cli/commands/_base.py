# coding: utf-8
import abc
import contextlib
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
import rich.progress
import rich.logging
from rich.text import Text

from ... import __version__, __name__ as __progname__

_T = typing.TypeVar("_T")


class InvalidArgument(ValueError):
    """An error to mark an invalid value was passed to a CLI flag.
    """

class CommandExit(Exception):
    """An error to request immediate exit from a function.
    """

    def __init__(self, code: int) -> None:
        self.code = code


class MofNWithPrecisionColumn(rich.progress.MofNCompleteColumn):

    def render(self, task: "rich.progress.Task") -> "Text":
        """Show completed/total."""
        completed = task.completed
        total = task.total if task.total is not None else "?"
        precision = task.fields.get("precision", "")
        total_width = len(f"{total:{precision}}")
        return Text(
            f"{completed:{total_width}{precision}}{self.separator}{total:{precision}}",
            style="progress.download",
        )


class Command(metaclass=abc.ABCMeta):
    """An abstract base class for ``gecco`` subcommands.
    """

    # -- Abstract methods ----------------------------------------------------

    summary: ClassVar[str] = NotImplemented

    @abc.abstractmethod
    def execute(self, ctx: contextlib.ExitStack) -> int:
        """Execute the command.

        Returns:
            `int`: The exit code for the command, with 0 on success, and any
            other number on error.

        """
        return NotImplemented

    @classmethod
    @abc.abstractmethod
    def doc(cls, fast: bool = False) -> str:
        """Get the help message for the command.
        """
        return NotImplemented

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

        self.progress = rich.progress.Progress(
            rich.progress.SpinnerColumn(finished_text="[green]:heavy_check_mark:[/]"),
            "[progress.description]{task.description}",
            rich.progress.BarColumn(bar_width=60),
            MofNWithPrecisionColumn(),
            "[progress.completed]{task.fields[unit]}",
            "[progress.percentage]{task.percentage:>3.0f}%",
            rich.progress.TimeElapsedColumn(),
            rich.progress.TimeRemainingColumn(),
            console=rich.console.Console(file=self.stream, soft_wrap=True),
            disable=self.quiet > 0,
            transient=True,
        )
        self.console = self.progress.console

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
        optional: bool = False,
        default: _T = None
    ) -> _T:
        _convert = (lambda x: x) if convert is None else convert
        _check = (lambda x: True) if check is None else check
        if self.args[name] is None:
            if optional:
                return default  # type: ignore
            self.error(f"Missing value for argument [purple]{name}[/]")
            raise InvalidArgument(self.args[name])
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
            return value  # type: ignore

    def _on_import_error(
        self,
        subcommand: str,
        e: ImportError
    ) -> None:
        import rich.traceback

        self.error(f"The [bold blue]{subcommand}[/] subcommand requires optional dependency [bold blue]{e.name}[/]")
        traceback = rich.traceback.Traceback.from_exception(type(e), e, e.__traceback__, extra_lines=0)
        self.console.print(traceback)


    # -- Logging methods -----------------------------------------------------

    def error(self, message: str, *args: Any, level: int = 0) -> None:
        if self.quiet <= 2 and level <= self.verbose:
            if self.verbose <= 1:
                self.console.print(
                    "[bold red]x[/]",
                    message,
                    *args,
                )
            else:
                self.console.print(
                    *self._logprefix(),
                    "[bold red]FAIL[/]",
                    message,
                    *args,
                )

    def info(self, verb: str, *args: Any, level: int = 1) -> None:
        if self.quiet == 0 and level <= self.verbose:
            if self.verbose <= 1:
                self.console.print(
                    "[bold blue]i[/]",
                    verb,
                    *args,
                )
            else:
                self.console.print(
                    *self._logprefix(),
                    f"[bold blue]INFO[/]",
                    verb,
                    *args,
                )

    def success(self, verb: str, *args: Any, level: int = 1) -> None:
        if self.quiet == 0 and level <= self.verbose:
            if self.verbose <= 1:
                self.console.print(
                    "[green]:heavy_check_mark:[/]",
                    verb,
                    *args,
                )
            else:
                self.console.print(
                    *self._logprefix(),
                    f"[bold green]  OK[/]",
                    verb,
                    *args,
                )

    def warn(self, verb: str, *args: Any, level: int = 0) -> None:
        if self.quiet <= 1 and level <= self.verbose:
            if self.verbose <= 1:
                self.console.print(
                    "[bold yellow]![/]",
                    verb,
                    *args,
                )
            else:
                self.console.print(
                    *self._logprefix(),
                    "[bold yellow]WARN[/]",
                    verb,
                    *args
                )

    def _logprefix(self) -> List[str]:
        return [
            f"[dim cyan]{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}[/]",
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
