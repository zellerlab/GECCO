import io
import os
import unittest
import tempfile
import shutil

from rich.console import Console

from gecco.cli import main

from ._base import TestCommand


class TestTrain(TestCommand, unittest.TestCase):
    name = "train"

    @property
    def folder(self):
        return os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_train_feature_type_domain(self):
        base = os.path.join(self.folder, "data", "mibig-2.0.proG2")
        clusters, features, genes = f"{base}.clusters.tsv", f"{base}.features.tsv", f"{base}.genes.tsv"

        argv = [
            "-vv", "train", "-f", features, "-c", clusters, "-o", self.tmpdir,
            "-g", genes, "--feature-type", "domain"
        ]
        with io.StringIO() as stream:
            retcode = main(argv, console=Console(file=stream))
            self.assertEqual(retcode, 0, stream.getvalue())

        files = os.listdir(self.tmpdir)
        self.assertIn("model.pkl", files)
        self.assertIn("model.pkl.md5", files)

    def test_train_feature_type_protein(self):
        base = os.path.join(self.folder, "data", "mibig-2.0.proG2")
        clusters, features, genes = f"{base}.clusters.tsv", f"{base}.features.tsv", f"{base}.genes.tsv"

        argv = [
            "-vv", "train", "-f", features, "-c", clusters, "-o", self.tmpdir,
            "-g", genes, "--feature-type", "protein"
        ]
        with io.StringIO() as stream:
            retcode = main(argv, console=Console(file=stream))
            self.assertEqual(retcode, 0, stream.getvalue())

        files = os.listdir(self.tmpdir)
        self.assertIn("model.pkl", files)
        self.assertIn("model.pkl.md5", files)
