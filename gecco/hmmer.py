"""Compatibility wrapper for HMMER binaries and output.
"""
import collections
import contextlib
import csv
import errno
import os
import subprocess
import tempfile
import typing
from typing import Dict, Optional, Iterable, List, Mapping, Type, Sequence

import pandas
from Bio import SeqIO

from ._base import BinaryRunner
from .model import Gene, Domain

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord
    from .data.hmms import Hmm

    _T = typing.TypeVar("_T", bound="DomainRow")


class DomainRow(typing.NamedTuple):
    """A single row in a domain table created by ``hmmsearch``.

    See Also:
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
    def from_line(cls: Type["_T"], row: str) -> "_T":
        """Extract a `DomainRow` from a single domain table line.
        """
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


class HMMER(BinaryRunner):
    """A wrapper for HMMER that scans a HMM library against protein sequences.
    """

    BINARY = "hmmsearch"

    def __init__(self, hmm: str, cpus: Optional[int] = None) -> None:
        """Prepare a new HMMER annotation handler with the given ``hmms``.

        Arguments:
            hmms (str): The path to the file containing the HMMs.
            cpus (int, optional): The number of CPUs to allocate for the
                ``hmmsearch`` command. Give ``None`` to use the default.

        """
        super().__init__()
        self.hmm = hmm
        self.cpus = cpus

    def run(self, genes: Sequence[Gene]) -> Sequence[Gene]:
        """Run HMMER on ``proteins`` and return found domains as a dataframe.

        Arguments:
            proteins (iterable of `~Bio.SeqRecord.SeqRecord`): The proteins to
                annotate with HMMER.
            prodigal (bool, optional): Whether or not the protein files were
                obtained with PRODIGAL, in which case the extraction of some
                features to the final dataframe will be a lot more accurate.
                Defaults to ``True``.

        """
        # collect genes and build an index of genes by protein id
        gene_index = { gene.id: gene for gene in genes }

        # create a temporary file to write the input and output to
        seqs_tmp = tempfile.NamedTemporaryFile(prefix="hmmer", suffix=".faa")
        doms_tmp = tempfile.NamedTemporaryFile(prefix="hmmer", suffix=".dom", mode="rt")

        # write protein sequences
        protein_records = (g.protein.to_record() for g in genes)
        SeqIO.write(protein_records, seqs_tmp.name, "fasta")

        # Prepare the command line arguments
        cmd = ["hmmsearch", "--noali", "--domtblout", doms_tmp.name]
        if self.cpus is not None:
            cmd.extend(["--cpu", str(self.cpus)])
        cmd.extend([self.hmm.path, seqs_tmp.name])

        # Run HMMER
        subprocess.run(cmd, stdout=subprocess.DEVNULL).check_returncode()

        # Read the domain table and update protein domains
        lines = filter(lambda line: not line.startswith("#"), doms_tmp)
        rows = map(DomainRow.from_line, lines)

        for row in rows:
            assert row.env_from <= row.env_to
            gene = gene_index[row.target_name]
            name = self.hmm.relabel(row.query_accession or row.query_name)
            domain = Domain(name, row.env_from, row.env_to, self.hmm.id, row.i_evalue)
            gene.protein.domains.append(domain)

        return genes
