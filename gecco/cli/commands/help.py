"""Implementation of the ``gecco run`` subcommand.
"""

import contextlib
import textwrap
from typing import Any, Dict, Mapping, List, Optional, TextIO

import rich.console
import rich.text

from .._utils import in_context, patch_showwarnings
from ._base import Command, CommandExit
from ._main import Main


class Help(Command):  # noqa: D101

    summary = "display the help message of another subcommand."

    @classmethod
    def doc(cls, fast=False):  # noqa: D102
        return f"""
        gecco help - {cls.summary}

        Usage:
            gecco help [<cmd>]

        Arguments:
            <cmd>                      a command to get the help message of.

        Parameters:
            -h, --help                 show the message for ``gecco`` or
                                       for a given subcommand.
        """

    def execute(self, ctx: contextlib.ExitStack) -> int:  # noqa: D102
        try:
            # check arguments and enter context
            self._check()
            ctx.enter_context(patch_showwarnings(self._showwarnings))
            # Get the subcommand class
            if self.args["<cmd>"] is not None:
                subcmd_cls = Main._get_subcommand_by_name(self.args["<cmd>"])
            else:
                subcmd_cls = None
            # Exit if no known command was found
            if self.args["<cmd>"] is not None and subcmd_cls is None:
                self.error("Unknown subcommand", repr(self.args["<cmd>"]))
                return 1
            # Render the help message
            doc = Main.doc() if subcmd_cls is None else subcmd_cls.doc()
            text = rich.text.Text(textwrap.dedent(doc).lstrip())
            console = rich.console.Console(file=self._stream, soft_wrap=True)
            console.print(text)
        except CommandExit as cexit:
            return cexit.code
        else:
            return 0
