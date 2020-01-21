import os
import errno
import subprocess

class ORFFinder(object):
    """Finds ORFs given a FASTA file and writes the results to out_dir"""

    def __init__(self, fasta, out_dir, method="prodigal",
            out_formats=None,
            other_args=None):
        self.fasta = fasta
        self.base = ".".join(os.path.basename(fasta).split(".")[:-1])
        self.out_dir = out_dir
        self.method = method
        self._check_method()

        self.out_formats =  ["genes", "coords"] if out_formats is None else out_formats
        self.other_args = [] if other_args is None else other_args

    def run(self):
        """Find ORFs"""
        cmd = self._make_commandline()
        log_out = os.path.join(self.out_dir, self.base + f".{self.method}.log")
        with open(log_out, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=out)
        return self.main_out

    def _check_method(self):
        """Checks wether chosen method is available. Raises error if not."""
        try:
            devnull = open(os.devnull)
            subprocess.run([self.method], stdout=devnull, stderr=devnull)
        except OSError as err:
            if err.errno == errno.ENOENT:
                raise RuntimeError(f"{self.method} does not seem to be installed. Please install it and re-run GECCO.") from None
            raise

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

        return cmd + ["-p", "meta", "-f", "sco"] + self.other_args
