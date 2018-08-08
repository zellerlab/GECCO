import os
import subprocess

class ORFFinder(object):
    """Finds ORFs given a FASTA file and writes the results to out_dir"""

    def __init__(self, fasta, out_dir, method="prodigal",
            out_formats=["genes", "coords"],
            other_args=["-p", "meta"]):
        self.fasta = fasta
        self.base = os.path.basename(self.fasta).split(".")[0]
        self.out_dir = out_dir
        self.method = method
        self._check_method()

        self.out_formats = out_formats
        self.other_args = other_args

    def run(self):
        """Find ORFs"""
        cmd = self._make_commandline()
        log_out = os.path.join(self.out_dir, self.base + f".{self.method}.log")
        subprocess.run(cmd,
            stdout = open(log_out, "wt"),
            stderr = open(log_out, "wt"))

        return self.main_out

    def _check_method(self):
        """Checks wether chosen method is available. Raises error if not."""
        try:
            devnull = open(os.devnull)
            subprocess.run([self.method], stdout=devnull, stderr=devnull)
        except OSError as err:
            if err.errno == os.errno.ENOENT:
                raise OSError(f"{self.method} does not seem to be installed. Please install it and re-run ORION.")

    def _make_commandline(self):
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

        return cmd + self.other_args
