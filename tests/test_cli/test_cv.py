import unittest
from gecco.cli.commands.cv import Cv
from ._base import TestCommand

class TestCv(TestCommand, unittest.TestCase):
    command_type = Cv
