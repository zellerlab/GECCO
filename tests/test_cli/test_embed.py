import unittest
from gecco.cli.commands.embed import Embed
from ._base import TestCommand

class TestEmbed(TestCommand, unittest.TestCase):
    command_type = Embed
