import csv
import errno
import os
import subprocess
import typing
from typing import Dict, Optional, List, Type

import pandas


_T = typing.TypeVar("_T", bound="DomainRow")


class DomainRow(typing.NamedTuple):
    """A single row in a domain table created by ``hmmsearch``.

    See also:
        The description of each field in page 48 of the `HMMER manual
        <http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf>`_.

    """

    target_name: str
    target_accession: Optional[str]
    target_length: int
    query_name: str
    query_accession: Optional[str]
    query_length: int
    evalue: float
    score: float
    bias: float
    domain_number: int
    domain_total: int
    c_evalue: float
    i_evalue: float
    domain_score: float
    domain_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    post: float
    description: Optional[str]

    @classmethod
    def from_line(cls: Type[_T], row: str) -> _T:
        line = list(filter(None, row.split(" ")))
        return cls(
            target_name=line[0],
            target_accession=None if line[1] == "-" else line[1],
            target_length=len(line[2]),
            query_name=line[3],
            query_accession=None if line[4] == "-" else line[4],
            query_length=len(line[5]),
            evalue=float(line[6]),
            score=float(line[7]),
            bias=float(line[8]),
            domain_number=int(line[9]),
            domain_total=int(line[10]),
            c_evalue=float(line[11]),
            i_evalue=float(line[12]),
            domain_score=float(line[13]),
            domain_bias=float(line[14]),
            hmm_from=int(line[15]),
            hmm_to=int(line[16]),
            ali_from=int(line[17]),
            ali_to=int(line[18]),
            env_from=int(line[19]),
            env_to=int(line[20]),
            post=float(line[21]),
            description=" ".join(line[22:]) if line[22:] else None,
        )


class HMMER(object):
    """A wrapper for HMMER that scans a HMM library against protein sequences.
    """

    def __init__(
        self,
        fasta: str,
        out_dir: str,
        hmms: str,
        prodigal: bool = True,
        cpus: Optional[int] = None,
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
            if e.errno == errno.ENOENT:
                raise OSError(
                    "HMMER does not seem to be installed. Please install it and re-run GECCO."
                )
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
            "domain_end",
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
                writer.writerow(
                    [sid, pid, start, end, strand, domain, l[12]] + typing.cast(List[object], l[17:19])
                )

    def _to_dataframe(self, dom_file: str) -> pandas.DataFrame:
        """Converts a HMMER domain table to a `pandas.DataFrame`.
        """
        rows = []
        with open(dom_file, "r") as f:
            for line in filter(lambda line: not line.startswith("#"), f):
                row = DomainRow.from_line(line)
                if self.prodigal:
                    sid = row.target_name[: row.target_name.rfind("_")]
                    pid = row.target_name
                    # extract additional metadata from the target description
                    description = typing.cast(str, row.description)
                    info = [x.strip() for x in description.split("#") if x]
                    start = int(info[0])
                    end = int(info[1])
                    strand = "+" if info[2] == "1" else "-"
                else:
                    sid = pid = row.target_name
                    start = self.protein_order[pid]
                    end = self.protein_order[pid]
                    strand = "unknown"
                domain = row.query_accession or row.query_name
                rows.append(
                    (
                        sid,
                        pid,
                        min(start, end),
                        max(start, end),
                        strand,
                        domain,
                        row.i_evalue,
                        1 - row.i_evalue,
                        min(row.env_from, row.env_to),
                        max(row.env_from, row.env_to),
                    )
                )
        return pandas.DataFrame(
            rows,
            columns=[
                "sequence_id",
                "protein_id",
                "start",
                "end",
                "strand",
                "domain",
                "i_Evalue",
                "rev_i_Evalue",
                "domain_start",
                "domain_end",
            ],
        )

    def _get_protein_order(self) -> Dict[str, int]:
        with open(self.fasta, "r") as f:
            pids = [line[1:].split()[0] for line in f if line.startswith(">")]
        return {pid: i for i, pid in enumerate(pids)}
