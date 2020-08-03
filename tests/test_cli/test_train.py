import unittest
from gecco.cli.commands.train import Train
from ._base import TestCommand

class TestTrain(TestCommand, unittest.TestCase):
    command_type = Train
