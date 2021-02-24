"""Implementation of the main ``gecco`` command.
"""

import sys
import textwrap
import typing
import warnings
from typing import Mapping, Optional, Type

import better_exceptions
import docopt
import operator
import pkg_resources
import rich.traceback

from ... import __version__
from .._utils import classproperty
from . import __name__ as __parent__
from ._base import Command


class Main(Command):
    """The *main* command launched before processing subcommands.
    """

    @classmethod
    def _get_subcommands(cls) -> Mapping[str, Type[Command]]:
        commands = {}
        for cmd in pkg_resources.iter_entry_points(__parent__):
            try:
                commands[cmd.name] = cmd.load()
            except pkg_resources.DistributionNotFound as err:
                pass
        return commands

    @classmethod
    def _get_subcommand(cls, name: str) -> Optional[Type[Command]]:
        return cls._get_subcommands().get(name)

    @classproperty
    def doc(cls) -> str:  # type: ignore
        commands = (
            "    {:12}{}".format(name, typing.cast(Command, cmd).summary)
            for name, cmd in sorted(
                cls._get_subcommands().items(), key=operator.itemgetter(0)
            )
        )
        return (
            textwrap.dedent(
                """
        gecco - Gene Cluster Prediction with Conditional Random Fields

        Usage:
            gecco [-v | -vv | -q | -qq] [options] (-h | --help) [<cmd>]
            gecco [-v | -vv | -q | -qq] [options] <cmd> [<args>...]

        Commands:
        {commands}

        Parameters:
            -h, --help                 show the message for ``gecco`` or
                                       for a given subcommand.
            -q, --quiet                silence any output other than errors
                                       (-qq silences everything).
            -v, --verbose              increase verbosity (-v is minimal,
                                       -vv is verbose, and -vvv shows
                                       debug information).
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
        retcode = self._check()
        if retcode is not None:
            return retcode

        # Get the subcommand class
        subcmd_name = self.args["<cmd>"]
        subcmd_cls = self._get_subcommand(self.args["<cmd>"])

        # Exit if no known command was found
        if subcmd_name is not None and subcmd_cls is None:
            self.error("Unknown subcommand", repr(subcmd_name))
            return 1

        # Print a help message if asked for
        if (
            self.args["--help"]
            or "-h" in self.args["<args>"]
            or "--help" in self.args["<args>"]
        ):
            subcmd = typing.cast(Type[Command], self._get_subcommand("help"))(
                argv=["help"] + [subcmd_name],
                stream=self._stream,
                logger=self.logger,
                options=self.args,
                config=self.config,
            )

        # Print version information
        elif self.args["--version"]:
            self.console.print("gecco", __version__)
            return 0

        # Initialize the command if is valid
        else:
            subcmd = typing.cast(Type[Command], subcmd_cls)(
                argv=[self.args["<cmd>"]] + self.args["<args>"],
                stream=self._stream,
                options=self.args,
                config=self.config,
            )
            subcmd.verbose = self.verbose
            subcmd.quiet = self.quiet

        # Run the app, elegantly catching any interrupts or exceptions
        try:
            exitcode = subcmd._check()
            if exitcode is None:
                exitcode = subcmd()
        except KeyboardInterrupt:
            self.error("interrupted")
            return 2
        except Exception as e:
            self.error(
                "An unexpected error occurred. Consider opening"
                " a new issue on the bug tracker"
                " (https://github.com/zellerlab/GECCO/issues/new) if"
                " it persists, including the traceback below:"
            )
            traceback = rich.traceback.Traceback.from_exception(type(e), e, e.__traceback__)
            self.console.print(traceback)
            # return errno if exception has any
            return typing.cast(int, getattr(e, "errno", 1))
        else:
            return exitcode
