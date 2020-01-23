# coding: utf-8
import textwrap
import typing

import better_exceptions
import docopt
import pkg_resources

from ... import __version__
from .._meta import classproperty, wrap_warnings
from ._base import Command
from . import __name__ as __parent__


class Main(Command):
    @classmethod
    def _get_subcommands(cls) -> typing.Mapping[str, Command]:
        return {
            cmd.name: cmd.load()
            for cmd in pkg_resources.iter_entry_points(__parent__)
        }

    @classmethod
    def _get_subcommand(cls, name: str) -> typing.Optional[Command]:
        try:
            return next(
                typing.cast(Command, ep.load())
                for ep in pkg_resources.iter_entry_points(__parent__)
                if ep.name == name
            )
        except StopIteration:
            return None

    @classproperty
    def doc(cls):
        commands = (
            "    {:27}{}".format(name, cmd.summary) for name, cmd in cls._get_subcommands().items()
        )
        return (
            textwrap.dedent(
                """
        gecco - gene cluster prediction with conditional random fields

        Usage:
            gecco [options] <cmd> [<args>...]
            gecco [options] (-h | --help) [<cmd>]

        Commands:
        {commands}

        Parameters:
            -h, --help                 show this message

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

        # Assert CLI arguments were parsed Successfully
        if isinstance(self.args, docopt.DocoptExit):
            print(self.args, file=self.stream)
            return 1

        # Get the subcommand class
        subcmd_cls = self._get_subcommand(self.args["<cmd>"])

        # Exit if no command was found
        if self.args["<cmd>"] is not None and subcmd_cls is None:
            self.logger.error("Unknown subcommand: {!r}", self.args["<cmd>"])
            return 1

        # Setup better exceptions if traceback is rendered
        if self.args["--traceback"]:
            better_exceptions.hook()

        # Print a help message if asked for
        if self.args["--help"]:
            doc = self.doc if subcmd_cls is None else subcmd_cls.doc
            print(textwrap.dedent(doc).lstrip(), file=self.stream)
            return 0

        # Initialize the command if is valid
        else:
            subcmd = wrap_warnings(self.logger)(
                wrap_warnings(self.logger)(
                    subcmd_cls(
                        argv=[self.args["<cmd>"]] + self.args["<args>"],
                        stream=self.stream,
                        logger=self.logger,
                        options=self.args,
                        config=self.config,
                    )
                )
            )

        # Run the app, elegantly catching any interrupts or exceptions
        try:
            exitcode = subcmd()
        except KeyboardInterrupt:
            self.logger.error("Interrupted")
            return 2
        except Exception as e:
            self.logger.critical("{}", e)
            if self.args["--traceback"]:
                raise
            return getattr(e, "errno", 1)  # return errno if exception has any
        else:
            return exitcode
