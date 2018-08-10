import os
import subprocess
import pandas as pd

class HMMER(object):
    """Searches a HMM library against protein sequences."""

    def __init__(self, fasta, out_dir, hmms):
        self.fasta = fasta
        self.out_dir = out_dir
        self.hmms = hmms
        self._check_hmmer()

        with open(fasta, "rt") as f:
            self.gene_order = [l.split()[0][1:] for l in f if l.startswith(">")]

    def run(self):
        """Runs HMMER. Converts output to tsv."""
        base = os.path.basename(self.fasta).split(".")[0]
        dom_out = os.path.join(self.out_dir, base + ".hmmer.dom")
        cmd =  ["hmmsearch", "--domtblout",
            dom_out, self.hmms, self.fasta]
        std_out = os.path.join(self.out_dir, base + ".hmmer.out")
        err_out = os.path.join(self.out_dir, base + ".hmmer.err")
        subprocess.run(cmd,
            stdout = open(std_out, "wt"),
            stderr = open(err_out, "wt"))

        tsv_out = os.path.join(self.out_dir, base + ".hmmer.tsv")
        self._to_tsv(dom_out, tsv_out)

        out_df = pd.read_csv(tsv_out, sep = "\t")
        out_df["protein_id"] = pd.Categorical(out_df["protein_id"], self.gene_order)
        out_df = out_df.sort_values(["protein_id", "domain_start"]).reset_index(drop=True)
        return out_df

    def _check_hmmer(self):
        """Checks wether hmmsearch is available. Raises error if not."""
        try:
            devnull = open(os.devnull)
            subprocess.run(["hmmsearch"], stdout=devnull, stderr=devnull)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise OSError("HMMER does not seem to be installed. Please install it and re-run ORION.")

    def _to_tsv(self, dom_file, out_file):
        """Converts HMMER --domtblout output to regular TSV"""
        header = ["protein_id", "pfam", "i_Evalue", "domain_start", "domain_end"]
        with open(dom_file, "rt") as f:
            with open(out_file, "wt") as fout:
                fout.write("\t".join(header) + "\n")
                for line in f:
                    if not line.startswith("#"):
                        line = line.split()
                        row = [line[0]] + [line[4]] + [line[12]] + line[17:19]
                        fout.write("\t".join(row) + "\n")
