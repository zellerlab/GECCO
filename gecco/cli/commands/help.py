import csv
import logging
import multiprocessing
import os
import pickle
import random
import textwrap
import typing

import numpy
import pandas
from Bio import SeqIO

from ._base import Command
from ._main import Main


class Help(Command):

    summary = "display the help message of another subcommand."
    doc = f"""
    gecco help - {summary}

    Usage:
        gecco help (-h | --help)
        gecco help [<cmd>]

    Arguments:
        <cmd>                      a command to get the help message of.

    Parameters:
        -h, --help                 show the message for ``gecco`` or
                                   for a given subcommand.
    """

    def __call__(self) -> int:
        # Get the subcommand class
        if self.args["<cmd>"] is not None:
            subcmd_cls = Main._get_subcommand(self.args["<cmd>"])
        else:
            subcmd_cls = None

        # Exit if no known command was found
        if self.args["<cmd>"] is not None and subcmd_cls is None:
            self.logger.error("Unknown subcommand: {!r}", self.args["<cmd>"])
            return 1

        # Render the help message
        doc = Main.doc if subcmd_cls is None else subcmd_cls.doc
        print(textwrap.dedent(doc).lstrip())
        return 0
