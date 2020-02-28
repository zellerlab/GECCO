import csv
import os
import subprocess
import typing

import pandas


class HMMER(object):
    """A wrapper for HMMER that scans a HMM library against protein sequences.
    """

    def __init__(
        self,
        fasta: str,
        out_dir: str,
        hmms: str,
        prodigal: bool = True,
        cpus: typing.Optional[int] = None,
    ) -> None:
        """Prepare a new `HMMER` annotation run.

        Arguments:
            fasta (str): The path to the file containing the input sequences.
            out_dir (str): The path to the directory in which to write output.
            hmms (str): The path to the file containing the HMMs.
            prodigal (bool, optional): Whether or not the protein files were
                obtained with PRODIGAL, in which case the extraction of some
                features to the final dataframe will be a lot more accurate.
                Defaults to ``True``.
            cpus (int, optional): The number of CPUs to allocate for the
                ``hmmsearch`` command. Give ``None`` to use the default.

        """
        self.fasta = fasta
        self.prodigal = prodigal
        if not self.prodigal:
            self.protein_order = self._get_protein_order()
        self.out_dir = out_dir
        self.hmms = hmms
        self.cpus = cpus
        self._check_hmmer()

    def run(self) -> pandas.DataFrame:
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

        # Extract the result as a dataframe
        return (
            self._to_dataframe(dom_out)
                .sort_values(["sequence_id", "start", "domain_start"])
                .reset_index(drop=True)
        )

    def _check_hmmer(self) -> None:
        """Checks wether hmmsearch is available. Raises error if not."""
        try:
            devnull = subprocess.DEVNULL
            subprocess.run(["hmmsearch"], stdout=devnull, stderr=devnull)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise OSError("HMMER does not seem to be installed. Please install it and re-run GECCO.")
            raise

    def _to_tsv(self, dom_file: str, out_file: str) -> None:
        """Converts HMMER --domtblout output to regular TSV"""
        header = [
            "sequence_id",
            "protein_id",
            "start",
            "end",
            "strand",
            "domain",
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
                    sid = "_".join(l[0].split("_")[:-1])
                    pid = l[0]
                    start = self.protein_order[pid]
                    end = self.protein_order[pid]
                    strand = "unknown"
                domain = l[4] or l[3]
                writer.writerow([sid, pid, start, end, strand, domain, l[12]] + l[17:19])

    def _to_dataframe(self, dom_file: str) -> pandas.DataFrame:
        """Converts a HMMER domain table to a `pandas.DataFrame`.
        """
        rows = []
        with open(dom_file, "r") as f:
            for line in filter(lambda line: not line.startswith("#"), f):
                l = list(filter(None, line.split(" ")))
                if self.prodigal:
                    sid = "_".join(l[0].split("_")[:-1])
                    pid = l[0]
                    start = int(l[23])
                    end = int(l[25])
                    strand = "+" if l[27] == "1" else "-"
                else:
                    sid = pid = l[0]
                    start = self.protein_order[pid]
                    end = self.protein_order[pid]
                    strand = "unknown"
                domain = l[3] if l[4] == "-" else l[4]
                domain_start, domain_end = int(l[17]), int(l[19])
                rows.append((
                    sid,
                    pid,
                    min(start, end),
                    max(start, end),
                    strand,
                    domain,
                    float(l[12]),
                    min(domain_start, domain_end),
                    max(domain_start, domain_end),
                ))
        return pandas.DataFrame(rows, columns=[
            "sequence_id", "protein_id", "start", "end", "strand",
            "domain", "i_Evalue", "domain_start", "domain_end",
        ])

    def _get_protein_order(self) -> typing.Dict[str, int]:
        with open(self.fasta, "r") as f:
            pids = [line[1:].split()[0] for line in f if line.startswith(">")]
        return {pid:i for i, pid in enumerate(pids)}
