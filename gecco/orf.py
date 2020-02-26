import errno
import os
import subprocess
import typing
from typing import Optional, List


class ORFFinder(object):
    """Finds ORFs given a FASTA file containing a genome.
    """

    def __init__(
            self,
            fasta: str,
            out_dir: str,
            method: str ="prodigal",
            out_formats: Optional[List[str]] = None,
            other_args: Optional[List[str]] = None
    ) -> None:
        """Create a new `ORFFinder` to launch PRODIGAL on ``fasta``.

        Arguments:
            fasta (str): The path to the FASTA file containing the input genome.
            out_dir (str): The path to the directory in which to write the
                results files.
            method (str): The method to use to extract the ORFs. Only
                *prodigal* is actually supported at the moment.
            out_formats (list, optional): A list of output formats to produce
                in supplement of the FASTA file containing the proteins.
                Defaults to ``["genes", "coords"]``.
            other_args (list, optional): A list of arbitrary arguments to pass
                to the invoked binary. Defaults to ``[]``,
        """

        self.fasta = fasta
        self.base, _ = os.path.splitext(os.path.basename(fasta))
        self.out_dir = out_dir
        self.method = method
        self._check_method()

        self.out_formats =  ["genes", "coords"] if out_formats is None else out_formats
        self.other_args: List[str] = [] if other_args is None else other_args

    def run(self) -> str:
        """Launch the ORF finder with the arguments passed at initialisation.

        Returns:
            `str`: the path to the file containing the translated ORFs that
            were discovered in the input file.
        """
        cmd = self._make_commandline()
        log_out = os.path.join(self.out_dir, self.base + f".{self.method}.log")
        with open(log_out, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=out).check_returncode()
        return self.main_out

    def _check_method(self) -> None:
        """Checks wether chosen method is available. Raises error if not."""
        try:
            devnull = subprocess.DEVNULL
            subprocess.run([self.method], stdout=devnull, stderr=devnull)
        except OSError as err:
            if err.errno == errno.ENOENT:
                raise RuntimeError(f"{self.method} does not seem to be installed. Please install it and re-run GECCO.") from None
            raise

    def _make_commandline(self) -> List[str]:
        """Makes commandline for subprocess"""
        cmd = [self.method, "-i", self.fasta]
        self.main_out = os.path.join(self.out_dir, self.base + ".proteins.faa")
        cmd += ["-a", self.main_out]

        if "genes" in self.out_formats:
            genes_out = os.path.join(self.out_dir, self.base + ".genes.fna")
            cmd += ["-d", genes_out]

        if "coords" in self.out_formats:
            coords_out = os.path.join(self.out_dir, self.base + ".coords.sco")
            cmd += ["-o", coords_out]

        return cmd + ["-p", "meta", "-f", "sco"] + self.other_args
