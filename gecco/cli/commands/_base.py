# coding: utf-8
import abc
import logging
import sys
import textwrap
import typing
from typing import Any, ClassVar, Optional, List, Mapping, Dict, TextIO

import coloredlogs
import docopt
import verboselogs

from ... import __version__, __name__ as __progname__
from .._utils import BraceAdapter


class Command(metaclass=abc.ABCMeta):
    """An abstract base class for ``gecco`` subcommands.
    """

    # -- Abstract methods ----------------------------------------------------

    doc: ClassVar[str] = NotImplemented
    summary: ClassVar[str] = NotImplemented

    @abc.abstractmethod
    def __call__(self) -> int:
        return NotImplemented  # type: ignore

    # -- Concrete methods ----------------------------------------------------

    _version: str = "{} {}".format(__progname__, __version__)
    _options_first: bool = False

    def __init__(
        self,
        argv: Optional[List[str]] = None,
        stream: Optional[TextIO] = None,
        logger: Optional[logging.Logger] = None,
        options: Optional[Mapping[str, Any]] = None,
        config: Optional[Dict[Any, Any]] = None,
    ) -> None:

        self.argv = argv
        self.stream: TextIO = stream or sys.stderr
        self.options = options or dict()
        self.pool = None
        self.config = config

        # Parse command line arguments
        try:
            self.args = docopt.docopt(
                textwrap.dedent(self.doc).lstrip(),
                help=False,
                argv=argv,
                version=self._version,
                options_first=self._options_first,
            )
            loglevel = self._get_log_level()
        except docopt.DocoptExit as de:
            self.args = de
            loglevel = None

        # Create a new colored logger if needed
        if logger is None:
            logger = verboselogs.VerboseLogger(__progname__)
            loglevel = (loglevel or "INFO").upper()
            coloredlogs.install(
                logger=logger,
                level=int(loglevel) if loglevel.isdigit() else loglevel,
                stream=self.stream,
            )

        # Use a loggin adapter to use new-style formatting
        self.logger = BraceAdapter(logger)

    def _check(self) -> Optional[int]:
        # Assert CLI arguments were parsed Successfully
        if isinstance(self.args, docopt.DocoptExit):
            print(self.args, file=self.stream)
            return 1
        # Display help if needed
        elif self.args["--help"]:
            print(textwrap.dedent(self.doc).lstrip())
            return 0
        else:
            return None

    def _get_log_level(self) -> Optional[str]:
        if self.args.get("--verbose"):
            return "VERBOSE" if self.args.get("--verbose") == 1 else "DEBUG"
        elif self.args.get("--quiet"):
            return "ERROR"
        else:
            return typing.cast(Optional[str], self.args.get("--log"))
