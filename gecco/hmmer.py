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

    def run(self):
        """Runs HMMER. Converts output to tsv."""
        base = ".".join(os.path.basename(self.fasta).split(".")[:-1])
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
        out_df = (out_df.sort_values(["sequence_id", "start", "domain_start"])
                        .reset_index(drop=True))
        return out_df

    def _check_hmmer(self):
        """Checks wether hmmsearch is available. Raises error if not."""
        try:
            devnull = open(os.devnull)
            subprocess.run(["hmmsearch"], stdout=devnull, stderr=devnull)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise OSError("HMMER does not seem to be installed. Please install it and re-run GECCO.")

    def _to_tsv(self, dom_file, out_file):
        """Converts HMMER --domtblout output to regular TSV"""
        header = ["sequence_id", "protein_id", "start", "end", "strand", "pfam",
            "i_Evalue", "domain_start", "domain_end"]
        with open(dom_file, "rt") as f:
            with open(out_file, "wt") as fout:
                fout.write("\t".join(header) + "\n")
                for l in f:
                    if not l.startswith("#"):
                        l = l.split()
                        sid = "_".join(l[0].split("_")[:-1])
                        start = min(l[23], l[25])
                        end = max(l[23], l[25])
                        strand = "+" if l[27] == "1" else "-"
                        row = [sid, l[0], start, end, strand, l[4], l[12]] + l[17:19]
                        fout.write("\t".join(row) + "\n")
