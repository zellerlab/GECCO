import abc
import datetime
import os
import socket
from typing import Any, List, Type, Optional, TextIO

import rich.progress
from rich.console import Console
from rich.progress import Progress, DownloadColumn

from .. import __name__ as PROGRAM
from ._utils import MofNCompleteUnitColumn


class ConsoleLogger(abc.ABC):
    console: Console

    def __init__(
        self,
        console: Console,
        *,
        quiet: int = 0,
        verbose: int = 0,
        hostname: Optional[str] = None,
        pid: Optional[int] = None,
        program: str = PROGRAM,
    ):
        self.console = console
        self.verbose = verbose
        self.quiet = quiet
        self._hostname = socket.gethostname() if hostname is None else hostname
        self._pid = os.getpid() if pid is None else pid
        self._program = program

    @abc.abstractmethod
    def error(self, message: str, *args: Any, level: int = 0) -> None:
        pass

    @abc.abstractmethod
    def info(self, verb: str, *args: Any, level: int = 1) -> None:
        pass

    @abc.abstractmethod
    def success(self, verb: str, *args: Any, level: int = 1) -> None:
        pass

    @abc.abstractmethod
    def warn(self, verb: str, *args: Any, level: int = 0) -> None:
        pass

    @abc.abstractmethod
    def progress(self, *args, **kwargs) -> Progress:
        pass


class FullConsoleLogger(ConsoleLogger):

    def _logprefix(self) -> List[str]:
        return [
            f"[bold dim cyan]{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}[/]",
            f"[bold dim purple]{self._hostname}[/]",
            f"[bold dim]{self._program}[[default dim]{self._pid}[/]][/]",
        ]

    def error(self, message: str, *args: Any, level: int = 0) -> None:
        if self.quiet <= 2 and level <= self.verbose:
            self.console.print(
                *self._logprefix(),
                "[bold red]FAIL[/]",
                message,
                *args,
            )

    def info(self, verb: str, *args: Any, level: int = 1) -> None:
        if self.quiet == 0 and level <= self.verbose:
            self.console.print(
                *self._logprefix(),
                f"[bold blue]INFO[/]",
                verb,
                *args,
            )

    def success(self, verb: str, *args: Any, level: int = 1) -> None:
        if self.quiet == 0 and level <= self.verbose:
            self.console.print(
                *self._logprefix(),
                f"[bold green]  OK[/]",
                verb,
                *args,
            )

    def warn(self, verb: str, *args: Any, level: int = 0) -> None:
        if self.quiet <= 1 and level <= self.verbose:
            self.console.print(*self._logprefix(), "[bold yellow]WARN[/]", verb, *args)

    def progress(self, download: bool = False) -> Progress:
        column = DownloadColumn() if download else MofNCompleteUnitColumn()
        return Progress(
            *self._logprefix(),
            f"[bold purple]WORK[/]",
            rich.progress.TextColumn("[progress.description]{task.description}"),
            rich.progress.BarColumn(),
            rich.progress.TaskProgressColumn(),
            rich.progress.TimeElapsedColumn(),
            # *args,
            # **kwargs,
            console=self.console,
            transient=True,
        )

class ConciseConsoleLogger(ConsoleLogger):

    def error(self, message: str, *args: Any, level: int = 0) -> None:
        if self.quiet <= 2 and level <= self.verbose:
            self.console.print(
                "[bold red]x[/]",
                message,
                *args,
            )

    def info(self, verb: str, *args: Any, level: int = 1) -> None:
        if self.quiet == 0 and level <= self.verbose:
            self.console.print(
                "[bold blue]i[/]",
                verb,
                *args,
            )

    def success(self, verb: str, *args: Any, level: int = 1) -> None:
        if self.quiet == 0 and level <= self.verbose:
            self.console.print(
                "[green]:heavy_check_mark:[/]",
                verb,
                *args,
            )

    def warn(self, verb: str, *args: Any, level: int = 0) -> None:
        if self.quiet <= 1 and level <= self.verbose:
            self.console.print(
                "[bold yellow]![/]",
                verb,
                *args,
            )

    def progress(self, download: bool = False) -> Progress:
        column = DownloadColumn() if download else MofNCompleteUnitColumn()
        return Progress(
            rich.progress.SpinnerColumn(finished_text="[green]:heavy_check_mark:[/]"),
            rich.progress.TextColumn("[progress.description]{task.description:<15}"),
            rich.progress.BarColumn(),
            rich.progress.TaskProgressColumn(),
            rich.progress.TimeElapsedColumn(),
            # *args,
            # **kwargs,
            column,
            console=self.console,
        )


def make_logger(
    console: Console,
    verbose: int,
    **kwargs,
):
    if verbose >= 2:
        return FullConsoleLogger(
            console=console,
            verbose=verbose,
            **kwargs,
        )
    else:
        return ConciseConsoleLogger(
            console=console,
            verbose=verbose,
            **kwargs,
        )


def showwarnings(
    logger: ConsoleLogger,
    message: str,
    category: Type[Warning],
    filename: str,
    lineno: int,
    file: Optional[TextIO] = None,
    line: Optional[str] = None,
    verbose: int = 0,
    quiet: int = 0,
) -> None:
    for line in filter(str.strip, str(message).splitlines()):
        logger.warn(line.strip())
