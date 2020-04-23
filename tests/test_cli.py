import io
import textwrap
import unittest
from unittest import mock

from gecco.cli import main
from gecco.cli.commands._main import Main
from gecco.cli.commands.annotate import Annotate
from gecco.cli.commands.cv import Cv
from gecco.cli.commands.embed import Embed
from gecco.cli.commands.help import Help
from gecco.cli.commands.run import Run
from gecco.cli.commands.train import Train
from gecco.cli.commands.tune import Tune


class TestCommand(object):

    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None
        cls.name = next(
            key
            for key, value in Main._get_subcommands().items()
            if value is cls.command_type
        )

    def test_help_redirection(self):
        # Check that, if a stream is given to `main`, it is used to print
        # everything, including help messages
        argv = [self.name, "--help"]
        stderr = io.StringIO()
        retcode = main(argv, stderr)
        self.assertEqual(retcode, 0)
        self.assertMultiLineEqual(
            stderr.getvalue().strip(),
            textwrap.dedent(self.command_type.doc).strip()
        )

    def test_help_stdout(self):
        # Check that when help is explicitly asked for, it gets printed to
        # stdout
        argv = [self.name, "--help"]
        with mock.patch("sys.stdout", new=io.StringIO()) as stdout:
            retcode = main(argv)
            self.assertEqual(retcode, 0)
            self.assertMultiLineEqual(
                stdout.getvalue().strip(),
                textwrap.dedent(self.command_type.doc).strip()
            )


class TestAnnotate(TestCommand, unittest.TestCase):
    command_type = Annotate

class TestCv(TestCommand, unittest.TestCase):
    command_type = Cv

class TestEmbed(TestCommand, unittest.TestCase):
    command_type = Embed

class TestHelp(TestCommand, unittest.TestCase):
    command_type = Help

class TestRun(TestCommand, unittest.TestCase):
    command_type = Run

class TestTrain(TestCommand, unittest.TestCase):
    command_type = Train

class TestTune(TestCommand, unittest.TestCase):
    command_type = Tune
