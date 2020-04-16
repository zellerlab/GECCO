# coding: utf-8
import abc
import logging
import sys
import textwrap
import typing

import coloredlogs
import docopt
import verboselogs

from ... import __version__, __name__ as __progname__
from .._meta import BraceAdapter


class Command(metaclass=abc.ABCMeta):
    """An abstract base class for ``gecco`` subcommands.
    """

    # -- Abstract methods ----------------------------------------------------

    doc: typing.ClassVar[str] = NotImplemented
    summary: typing.ClassVar[str] = NotImplemented

    @abc.abstractmethod
    def __call__(self) -> int:
        return NotImplemented  # type: ignore

    # -- Concrete methods ----------------------------------------------------

    _version: str = "{} {}".format(__progname__, __version__)
    _options_first: bool = False

    def __init__(
        self,
        argv: typing.Optional[typing.List[str]] = None,
        stream: typing.Optional[typing.TextIO] = None,
        logger: typing.Optional[logging.Logger] = None,
        options: typing.Optional[typing.Mapping[str, typing.Any]] = None,
        config: typing.Optional[typing.Dict[typing.Any, typing.Any]] = None,
    ) -> None:

        self.argv = argv
        self.stream: typing.TextIO = stream or sys.stderr
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
                level=int(loglevel) if loglevel.isdigit() else loglevel,
                stream=stream,
                logger=logger,
            )

        # Use a loggin adapter to use new-style formatting
        self.logger = BraceAdapter(logger)

    def _check(self) -> typing.Optional[int]:
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

    def _get_log_level(self) -> typing.Optional[str]:
        if self.args.get("--verbose"):
            return "VERBOSE" if self.args.get("--verbose") == 1 else "DEBUG"
        elif self.args.get("--quiet"):
            return "ERROR"
        else:
            return self.args.get("--log")
