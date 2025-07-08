import unittest
from ._base import TestCommand

class TestCv(TestCommand, unittest.TestCase):
    name = "cv"
