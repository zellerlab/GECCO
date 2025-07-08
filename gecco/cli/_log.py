import abc
import datetime
import os
import socket
from typing import Any, List, Type, Optional, TextIO

from .. import __name__ as PROGRAM
from rich.console import Console


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


class FullConsoleLogger(ConsoleLogger):

    def _logprefix(self) -> List[str]:
        return [
            f"[dim cyan]{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}[/]",
            f"[dim purple]{self._hostname}[/]",
            f"[dim]{self._program}[[default dim]{self._pid}[/]][/]",
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
            self.console.print(
                *self._logprefix(), 
                "[bold yellow]WARN[/]", 
                verb, 
                *args
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
