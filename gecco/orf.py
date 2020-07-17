"""Generic protocol for ORF detection in DNA sequences.
"""

import abc
import io
import os
import subprocess
import tempfile
import typing
from typing import Iterable, Iterator, List, Optional

import Bio.Alphabet
import Bio.SeqIO
import pyrodigal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ._base import BinaryRunner
from .model import Gene, Protein, Strand

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


class ORFFinder(metaclass=abc.ABCMeta):
    """An abstract base class to provide a generic ORF finder.
    """

    @abc.abstractmethod
    def find_genes(self, dna: SeqRecord) -> Iterable[Gene]:  # type: ignore
        """Find all genes from a DNA sequence.
        """
        return NotImplemented  # type: ignore


class PyrodigalFinder(ORFFinder):
    """An `ORFFinder` that uses the Pyrodigal bindings to Prodigal.

    Prodigal is a fast and reliable protein-coding gene prediction
    for prokaryotic genomes, with support for draft genomes and metagenomes.

    See Also:
        `Doug Hyatt, Gwo-Liang Chen, Philip F. LoCascio, Miriam L. Land,
        Frank W. Larimer and Loren J. Hauser.
        "Prodigal: Prokaryotic Gene Recognition and Translation Initiation
        Site Identification", BMC Bioinformatics 11 (8 March 2010), p119
        <https://doi.org/10.1186/1471-2105-11-119>`_

    """

    def __init__(self, metagenome: bool = True) -> None:
        """Create a new `PyrodigalFinder` instance.

        Arguments:
            metagenome (bool): Whether or not to run PRODIGAL in metagenome
                mode, defaults to `True`.

        """
        super().__init__()
        self.metagenome = metagenome
        self.pyrodigal = pyrodigal.Pyrodigal(meta=metagenome)

    def find_genes(self, dna: SeqRecord) -> Iterator[Gene]:  # noqa: D102
        # find all ORFs in the given DNA sequence
        orfs = self.pyrodigal.find_genes(str(dna.seq))
        for j, orf in enumerate(orfs):
            # wrap the protein into a Protein object
            prot_seq = Seq(orf.translate(), Bio.Alphabet.generic_protein)
            protein = Protein(id=f"{dna.id}_{j+1}", seq=prot_seq)
            # wrap the gene into a Gene
            yield Gene(
                source=dna,
                start=min(orf.begin, orf.end),
                end=max(orf.begin, orf.end),
                strand=Strand.Coding if orf.strand == 1 else Strand.Reverse,
                protein=protein,
            )
