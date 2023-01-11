"""Implementation of the main ``gecco`` command.
"""

import contextlib
import signal
import sys
import textwrap
import typing
import warnings
from typing import Mapping, Optional, Type, List

import docopt
import operator

from ... import __version__
from ..._meta import classproperty
from .._utils import in_context, patch_showwarnings
from . import __name__ as __parent__
from ._base import Command, CommandExit, InvalidArgument

if sys.version_info < (3, 10):
    import importlib_metadata  # type: ignore
else:
    import importlib.metadata as importlib_metadata


class Main(Command):
    """The *main* command launched before processing subcommands.
    """

    _entry_points_cache: Optional[List["importlib_metadata.EntryPoint"]] = None

    @classmethod
    def _entry_points(cls) -> List["importlib_metadata.EntryPoint"]:
        if cls._entry_points_cache is None:
            try:
                entry_points = importlib_metadata.entry_points(group=__parent__)
                cls._entry_points_cache = list(entry_points)
            except KeyError:
                cls._entry_points_cache = []
        return cls._entry_points_cache

    @classmethod
    def _get_subcommand_names(cls) -> List[str]:
        return [cmd.name for cmd in cls._entry_points()]

    @classmethod
    def _get_subcommands(cls) -> Mapping[str, Type[Command]]:
        commands = {}
        for cmd in cls._entry_points():
            try:
                commands[cmd.name] = cmd.load()
            except Exception:
                pass
        return commands

    @classmethod
    def _get_subcommand_by_name(cls, name: str) -> Optional[Type[Command]]:
        for cmd in cls._entry_points():
            if cmd.name == name:
                return cmd.load()  # type: ignore
        return None

    # --

    @classmethod
    def doc(cls, fast: bool = False) -> str:  # noqa: D102
        if fast:
            commands = (f"    {cmd}" for cmd in cls._get_subcommand_names())
        else:
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
            gecco [-v | -vv | -q | -qq] [--no-progress-bar] <cmd> [<args>...]
            gecco --version
            gecco --help [<cmd>]

        Commands:
        {commands}

        Parameters:
            -h, --help                 show the message for ``gecco`` or
                                       for a given subcommand.
            -q, --quiet                silence any output other than errors
                                       (-qq silences everything).
            -v, --verbose              increase verbosity (-v is verbose,
                                       -vv is very verbose and makes the
                                       output more suitable for logging).
            --no-progress-bar          disable progress bars from the output,
                                       independently of verbosity.
            -V, --version              show the program version and exit.

        """
            )
            .lstrip()
            .format(commands="\n".join(commands))
        )

    _options_first = True

    # --

    def execute(self, ctx: contextlib.ExitStack) -> int:
        # Run the app, elegantly catching any interrupts or exceptions
        try:
            # check arguments and enter context
            self._check()
            ctx.enter_context(patch_showwarnings(self._showwarnings)) # type: ignore

            # Get the subcommand class
            subcmd_name = self.args["<cmd>"]
            try:
                subcmd_cls = self._get_subcommand_by_name(subcmd_name)
            except ImportError as err:
                self._on_import_error(subcmd_name, err)
                return 1

            # exit if no known command was found
            if subcmd_name is not None and subcmd_cls is None:
                self.error("Unknown subcommand", repr(subcmd_name))
                return 1
            # if a help message was required, delegate to the `gecco help` command
            if (
                self.args["--help"]
                or "-h" in self.args["<args>"]
                or "--help" in self.args["<args>"]
            ):
                subcmd = typing.cast(Type[Command], self._get_subcommand_by_name("help"))(
                    argv=["help"] + [subcmd_name],
                    stream=self._stream,
                    options=self.args,
                    config=self.config,
                )
            # print version information if `--version` in flags
            elif self.args["--version"]:
                self.console.print("gecco", __version__)
                return 0
            # initialize the command if is valid
            else:
                subcmd = typing.cast(Type[Command], subcmd_cls)(
                    argv=[self.args["<cmd>"]] + self.args["<args>"],
                    stream=self._stream,
                    options=self.args,
                    config=self.config,
                )
                subcmd.verbose = self.verbose
                subcmd.quiet = self.quiet
                subcmd.progress.disable = self.args["--no-progress-bar"]
            # run the subcommand
            return subcmd.execute(ctx)
        except CommandExit as sysexit:
            return sysexit.code
        except KeyboardInterrupt:
            self.error("Interrupted")
            return -signal.SIGINT
        except ImportError as err:
            self._on_import_error(subcmd_name, err)
            return 1
        except Exception as e:
            import rich.traceback

            self.error(
                "An unexpected error occurred. Consider opening"
                " a new issue on the bug tracker"
                " ( https://github.com/zellerlab/GECCO/issues/new ) if"
                " it persists, including the traceback below:"
            )
            traceback = rich.traceback.Traceback.from_exception(type(e), e, e.__traceback__)
            self.console.print(traceback)
            # return errno if exception has any
            return typing.cast(int, getattr(e, "errno", 1))
