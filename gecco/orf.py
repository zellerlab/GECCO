import abc
import errno
import os
import subprocess
import tempfile
import typing
from typing import Optional, List

import Bio.SeqIO

from ._base import BinaryRunner

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


class ORFFinder(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def find_proteins(
        self,
        sequences: List["SeqRecord"],
    ) -> List["SeqRecord"]:
        return NotImplemented  # tupe: ignore


class ProdigalFinder(BinaryRunner, ORFFinder):

    BINARY = "prodigal"

    def __init__(self, metagenome=True):
        super().__init__()
        self.metagenome = metagenome

    def find_proteins(
        self,
        sequences: List["SeqRecord"],
    ) -> List["SeqRecord"]:
        _, seqs_path = tempfile.mkstemp(prefix=self.BINARY, suffix=".fna")
        _, prot_path = tempfile.mkstemp(prefix=self.BINARY, suffix=".faa")

        #
        try:
            # write a FASTA file containing the sequences in a temporary location
            Bio.SeqIO.write(sequences, seqs_path, "fasta")
            # build the command line
            cmd: List[str] = [self.BINARY, "-q", "-i", seqs_path, "-a", prot_path]
            if self.metagenome:
                cmd.extend(["-p", "meta"])

            subprocess.run(cmd, stdout=subprocess.DEVNULL).check_returncode()
            return list(Bio.SeqIO.parse(prot_path, "fasta"))
        finally:
            os.remove(seqs_path)
            os.remove(prot_path)
