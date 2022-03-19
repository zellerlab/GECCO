import io
import os
import unittest
import tempfile
import shutil

from gecco.cli import main
from gecco.cli.commands.train import Train

from ._base import TestCommand


class TestTrain(TestCommand, unittest.TestCase):
    command_type = Train

    @property
    def folder(self):
        return os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_train_embedding(self):
        base = os.path.join(self.folder, "data", "mibig-2.0.proG2")
        clusters, features, genes = f"{base}.clusters.tsv", f"{base}.features.tsv", f"{base}.gff"

        argv = ["-vv", "train", "-f", features, "-c", clusters, "-o", self.tmpdir, "-g", genes]
        with io.StringIO() as stream:
            retcode = main(argv, stream=stream)
            self.assertEqual(retcode, 0, stream.getvalue())

        files = os.listdir(self.tmpdir)
        self.assertIn("model.pkl", files)
        self.assertIn("model.pkl.md5", files)
