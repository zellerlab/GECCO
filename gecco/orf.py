import abc
import io
import os
import subprocess
import tempfile
import typing
from typing import Iterable, List, Optional

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
        return NotImplemented  # type: ignore


class ProdigalFinder(BinaryRunner, ORFFinder):

    BINARY = "prodigal"

    def __init__(self, metagenome: bool = True) -> None:
        super().__init__()
        self.metagenome = metagenome

    def find_proteins(
        self,
        sequences: Iterable["SeqRecord"],
    ) -> List["SeqRecord"]:
        #
        with tempfile.NamedTemporaryFile("w+", prefix=self.BINARY, suffix=".faa") as tmp:
            # write a FASTA buffer to pass as PRODIGAL input
            buffer = io.TextIOWrapper(io.BytesIO())
            Bio.SeqIO.write(sequences, buffer, "fasta")
            # build the command line
            cmd: List[str] = [self.BINARY, "-q", "-a", tmp.name]
            if self.metagenome:
                cmd.extend(["-p", "meta"])
            # run the program
            completed = subprocess.run(
                cmd,
                input=buffer.detach().getbuffer(),
                stdout=subprocess.DEVNULL
            )
            completed.check_returncode()
            return list(Bio.SeqIO.parse(tmp, "fasta"))
