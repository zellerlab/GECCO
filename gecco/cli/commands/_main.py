"""Implementation of the main ``gecco`` command.
"""

import textwrap
import typing
from typing import Mapping, Optional, Type

import better_exceptions
import docopt
import operator
import pkg_resources

from ... import __version__
from .._utils import classproperty, wrap_warnings
from . import __name__ as __parent__
from ._base import Command


class Main(Command):
    """The *main* command launched before processing subcommands.
    """

    @classmethod
    def _get_subcommands(cls) -> Mapping[str, Type[Command]]:
        return {
            cmd.name: cmd.load()
            for cmd in pkg_resources.iter_entry_points(__parent__)
        }

    @classmethod
    def _get_subcommand(cls, name: str) -> Optional[Type[Command]]:
        try:
            return next(
                typing.cast(Type[Command], ep.load())
                for ep in pkg_resources.iter_entry_points(__parent__)
                if ep.name == name
            )
        except StopIteration:
            return None

    @classproperty
    def doc(cls) -> str: # type: ignore
        commands = (
            "    {:12}{}".format(name, typing.cast(Command, cmd).summary)
            for name, cmd in sorted(
                cls._get_subcommands().items(),
                key=operator.itemgetter(0)
            )
        )
        return (
            textwrap.dedent(
                """
        gecco - Gene Cluster Prediction with Conditional Random Fields

        Usage:
            gecco [-v | -vv | -q | -l <level>] [options] (-h | --help) [<cmd>]
            gecco [-v | -vv | -q | -l <level>] [options] <cmd> [<args>...]

        Commands:
        {commands}

        Parameters:
            -h, --help                 show the message for ``gecco`` or
                                       for a given subcommand.
            -q, --quiet                silence any output other than errors
                                       (corresponds to the log level ERROR).
            -v, --verbose              control the verbosity of the output
                                       (corresponds to the log level DEBUG).
            -V, --version              show the program version and exit.

        Parameters - Debug:
            --traceback                display full traceback on error.
            -l <level>, --log <level>  the level of log message to display.
                                       [available: DEBUG, INFO, WARNING, ERROR]
        """
            )
            .lstrip()
            .format(commands="\n".join(commands))
        )

    _options_first = True

    def __call__(self) -> int:
        # Assert CLI arguments were parsed successfully
        if isinstance(self.args, docopt.DocoptExit):
            print(self.args, file=self.stream)
            return 1

        # Get the subcommand class
        subcmd_cls = self._get_subcommand(self.args["<cmd>"])

        # Exit if no known command was found
        if self.args["<cmd>"] is not None and subcmd_cls is None:
            self.logger.error("Unknown subcommand: {!r}", self.args["<cmd>"])
            return 1

        # Setup better exceptions if traceback is rendered
        if self.args["--traceback"]:
            better_exceptions.hook()

        # Print a help message if asked for
        if self.args["--help"] or "-h" in self.args["<args>"] or "--help" in self.args["<args>"]:
            subcmd = typing.cast(Type[Command], self._get_subcommand("help"))(
                argv=["help"] + [self.args["<cmd>"]],
                stream=self._stream,
                logger=self.logger,
                options=self.args,
                config=self.config,
            )

        # Print version information
        elif self.args["--version"]:
            print("gecco", __version__)
            return 0

        # Initialize the command if is valid
        else:
            subcmd = wrap_warnings(self.logger)( # type: ignore
                typing.cast(Type[Command], subcmd_cls)(
                    argv=[self.args["<cmd>"]] + self.args["<args>"],
                    stream=self._stream,
                    logger=self.logger,
                    options=self.args,
                    config=self.config,
                )
            )

        # Run the app, elegantly catching any interrupts or exceptions
        try:
            exitcode = subcmd._check()
            if exitcode is None:
                exitcode = subcmd()
        except KeyboardInterrupt:
            self.logger.error("Interrupted")
            return 2
        except Exception as e:
            self.logger.critical("{}", e)
            if self.args["--traceback"]:
                raise
            # return errno if exception has any
            return typing.cast(int, getattr(e, "errno", 1))
        else:
            return exitcode
