import csv
import os
import subprocess
import typing

import pandas as pd


class HMMER(object):
    """Searches a HMM library against protein sequences.
    """

    def __init__(
            self,
            fasta: str,
            out_dir: str,
            hmms: str,
            prodigal: bool = True,
            cpus: typing.Optional[int] = None,
    ) -> None:
        self.fasta = fasta
        self.prodigal = prodigal
        if not self.prodigal:
            self.protein_order = self._get_protein_order()
        self.out_dir = out_dir
        self.hmms = hmms
        self.cpus = cpus
        self._check_hmmer()

    def run(self):
        """Runs HMMER and returns the output as a data frame.
        """
        base, _ = os.path.splitext(os.path.basename(self.fasta))
        dom_out = os.path.join(self.out_dir, f"{base}.hmmer.dom")
        stdout = os.path.join(self.out_dir, f"{base}.hmmer.out")
        stderr = os.path.join(self.out_dir, f"{base}.hmmer.err")

        # Prepare the command line arguments
        cmd = ["hmmsearch", "-o", stdout, "--domtblout", dom_out]
        if self.cpus is not None:
            cmd.extend(["--cpu", str(self.cpus)])
        cmd.extend([self.hmms, self.fasta])

        # Run HMMER
        with open(stderr, "w") as err:
            subprocess.run(cmd, stderr=err).check_returncode()

        # Convert to TSV
        tsv_out = os.path.join(self.out_dir, f"{base}.hmmer.tsv")
        self._to_tsv(dom_out, tsv_out)

        # Sort table properly
        out_df = pd.read_csv(tsv_out, sep = "\t")
        out_df = (out_df.sort_values(["sequence_id", "start", "domain_start"])
                        .reset_index(drop=True))
        return out_df

    def _check_hmmer(self):
        """Checks wether hmmsearch is available. Raises error if not."""
        try:
            devnull = subprocess.DEVNULL
            subprocess.run(["hmmsearch"], stdout=devnull, stderr=devnull)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise OSError("HMMER does not seem to be installed. Please install it and re-run GECCO.")

    def _to_tsv(self, dom_file: str, out_file: str) -> None:
        """Converts HMMER --domtblout output to regular TSV"""
        header = [
            "sequence_id",
            "protein_id",
            "start",
            "end",
            "strand",
            "pfam",
            "i_Evalue",
            "domain_start",
            "domain_end"
        ]
        with open(dom_file, "r") as f, open(out_file, "w") as fout:
            writer = csv.writer(fout, dialect="excel-tab")
            writer.writerow(header)

            for line in filter(lambda line: not line.startswith("#"), f):
                l = line.split()
                if self.prodigal:
                    sid = "_".join(l[0].split("_")[:-1])
                    pid = l[0]
                    start = min(int(l[23]), int(l[25]))
                    end = max(int(l[23]), int(l[25]))
                    strand = "+" if l[27] == "1" else "-"
                else:
                    sid = "NA"
                    pid = l[0]
                    start = self.protein_order[pid]
                    end = self.protein_order[pid]
                    strand = "NA"
                pfam = l[4] or l[3]
                writer.writerow([sid, pid, start, end, strand, pfam, l[12]] + l[17:19])

    def _get_protein_order(self):
        with open(self.fasta, "r") as f:
            pids = [line[1:].split()[0] for line in f if line.startswith(">")]
        return dict(enumerate(pids))
