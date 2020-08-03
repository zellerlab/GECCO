import unittest
from gecco.cli.commands.help import Help
from ._base import TestCommand

class TestHelp(TestCommand, unittest.TestCase):
    command_type = Help
