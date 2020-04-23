"""Generic protocol for ORF detection in DNA sequences.
"""

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
    """An abstract base class to provide a generic ORF finder.
    """

    @abc.abstractmethod
    def find_proteins(self, sequences: List["SeqRecord"],) -> List["SeqRecord"]:
        """Find all proteins from a list of DNA sequences.
        """
        return NotImplemented  # type: ignore


class ProdigalFinder(BinaryRunner, ORFFinder):
    """An `ORFFinder` that wraps the PRODIGAL binary.

    PRODIGAL is a fast and reliable protein-coding gene prediction for
    prokaryotic genomes, with support for draft genomes and metagenomes.

    See Also:
        .. [PMC2848648] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/

    """

    BINARY = "prodigal"

    def __init__(self, metagenome: bool = True) -> None:
        """Create a new `ProdigalFinder` instance.

        Arguments:
            metagenome (bool): Whether or not to run PRODIGAL in metagenome
                mode, defaults to `True`.

        """
        super().__init__()
        self.metagenome = metagenome

    def find_proteins(
        self, sequences: Iterable["SeqRecord"],
    ) -> List["SeqRecord"]:  # noqa: D102
        with tempfile.NamedTemporaryFile(
            "w+", prefix=self.BINARY, suffix=".faa"
        ) as tmp:
            # write a FASTA buffer to pass as PRODIGAL input
            buffer = io.TextIOWrapper(io.BytesIO())
            Bio.SeqIO.write(sequences, buffer, "fasta")
            # build the command line
            cmd: List[str] = [self.BINARY, "-q", "-a", tmp.name]
            if self.metagenome:
                cmd.extend(["-p", "meta"])
            # run the program
            completed = subprocess.run(
                cmd, input=buffer.detach().getbuffer(), stdout=subprocess.DEVNULL
            )
            completed.check_returncode()
            return list(Bio.SeqIO.parse(tmp, "fasta"))
