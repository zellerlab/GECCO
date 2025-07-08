import io
import textwrap
from unittest import mock

from rich.console import Console

from gecco.cli import main


class TestCommand(object):
    name = None

    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None

    def test_help_redirection(self):
        # Check that, if a stream is given to `main`, it is used to print
        # everything, including help messages
        argv = [self.name, "--help"]
        stderr = io.StringIO()
        console = Console(file=stderr, force_terminal=False)
        retcode = main(argv, console)
        self.assertEqual(retcode, 0)
        self.assertTrue(
            stderr.getvalue().strip().startswith(f"Usage: gecco {self.name}")
            # textwrap.dedent(self.command_type.doc()).strip()
        )

    def test_help_stdout(self):
        # Check that when help is explicitly asked for, it gets printed to
        # stdout
        argv = [self.name, "--help"]
        with mock.patch("sys.stdout", new=io.StringIO()) as stdout:
            retcode = main(argv)
            self.assertEqual(retcode, 0)
            self.assertTrue(
                stdout.getvalue().strip().startswith(f"Usage: gecco {self.name}")
            )
