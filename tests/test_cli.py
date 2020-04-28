import contextlib
import io
import os
import shutil
import tempfile
import textwrap
import unittest
from unittest import mock

import pandas

import gecco.cli.commands.run
from gecco.cli import main
from gecco.cli.commands._main import Main
from gecco.cli.commands.annotate import Annotate
from gecco.cli.commands.cv import Cv
from gecco.cli.commands.embed import Embed
from gecco.cli.commands.help import Help
from gecco.cli.commands.run import Run
from gecco.cli.commands.train import Train
from gecco.cli.commands.tune import Tune


DATADIR = os.path.realpath(os.path.join(__file__, "..", "data"))


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

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_fasta_genome(self):
        sequence = os.path.join(DATADIR, "BGC0001866.fna")
        domains = pandas.read_table(os.path.join(DATADIR, "BGC0001866.features.tsv"))

        # we mock time consuming operations (type prediction and HMM annotation)
        # with precomputed or fake results
        _run = mock.Mock(return_value=domains)
        _fit_predict = lambda self,tc,nc,y: [("Polyketide", 1) for _ in nc]
        with contextlib.ExitStack() as stack:
            stack.enter_context(
                mock.patch.object(gecco.cli.commands.run.HMMER, "run", new=_run)
            )
            stack.enter_context(
                mock.patch.object(gecco.cli.commands.run.ClusterKNN, "fit_predict", new=_fit_predict)
            )
            argv = ["-vv", "--traceback", "run", "--genome", sequence, "--output", self.tmpdir]
            main(argv, stream=io.StringIO())

        # make sure we have generated the files we want
        # and that we found one cluster
        output = os.listdir(self.tmpdir)
        self.assertIn("BGC0001866.features.tsv", output)
        self.assertIn("BGC0001866.clusters.tsv", output)
        clusters = pandas.read_table(os.path.join(self.tmpdir, "BGC0001866.clusters.tsv"))
        self.assertEqual(len(clusters), 1)


class TestTrain(TestCommand, unittest.TestCase):
    command_type = Train


class TestTune(TestCommand, unittest.TestCase):
    command_type = Tune
