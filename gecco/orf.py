"""Generic protocol for ORF detection in DNA sequences.
"""

import abc
import io
import itertools
import os
import queue
import tempfile
import typing
import multiprocessing.pool
from multiprocessing.sharedctypes import Value
from typing import Callable, Iterable, Iterator, List, Optional

import Bio.SeqIO
import pyrodigal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .model import Gene, Protein, Strand

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


__all__ = ["ORFFinder", "PyrodigalFinder"]


class ORFFinder(metaclass=abc.ABCMeta):
    """An abstract base class to provide a generic ORF finder.
    """

    @abc.abstractmethod
    def find_genes(self, records: Iterable[SeqRecord]) -> Iterable[Gene]:  # type: ignore
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

    def __init__(self, metagenome: bool = True, mask: bool = False, cpus: int = 0) -> None:
        """Create a new `PyrodigalFinder` instance.

        Arguments:
            metagenome (bool): Whether or not to run Prodigal in metagenome
                mode, defaults to `True`.
            mask (bool): Whether or not to mask genes running across regions
                containing unknown nucleotides, defaults to `False`.
            cpus (int): The number of threads to use to run Pyrodigal in
                parallel. Pass ``0`` to use the number of CPUs on the machine.

        """
        super().__init__()
        self.metagenome = metagenome
        self.mask = mask
        self.cpus = cpus
        self.orf_finder =  pyrodigal.OrfFinder(meta=metagenome, mask=mask)

    def _train(self, records):
        sequences = []
        for i, seq in enumerate(parse(args.i)):
            if i > 0:
                sequences.append("TTAATTAATTAA")
            sequences.append(seq.seq)
        if len(sequences) > 1:
            sequences.append("TTAATTAATTAA")
        self.orf_finder.train("".join(sequences))

    def _process_record(self, record: SeqRecord) -> pyrodigal.Genes:
        return record, self.orf_finder.find_genes(str(record.seq))

    def find_genes(
        self,
        records: Iterable[SeqRecord],
        progress: Optional[Callable[[SeqRecord, int], None]] = None,
        *,
        pool_factory = multiprocessing.pool.ThreadPool,
    ) -> Iterator[Gene]:
        """Find all genes contained in a sequence of DNA records.

        Arguments:
            records (iterable of `~Bio.SeqRecord.SeqRecord`): An iterable of
                DNA records in which to find genes.
            progress (callable, optional): A progress callback of signature
                ``progress(record, total)`` that will be called everytime a
                record has been processed successfully, with ``record`` being
                the `~Bio.SeqRecord.SeqRecord` instance, and ``total`` being
                the total number of records to process.

        Keyword Arguments:
            pool_factory (`type`): The callable for creating pools, defaults
                to the `multiprocessing.pool.ThreadPool` class, but
                `multiprocessing.pool.Pool` is also supported.

        Yields:
            `~gecco.model.Gene`: An iterator over all the genes found in
            the given records.

        """
        # detect the number of CPUs
        _cpus = self.cpus if self.cpus > 0 else os.cpu_count()

        # train first if needed
        if not self.metagenome:
            records = list(records)
            self._train(records)

        # run in parallel using a pool
        with pool_factory(_cpus) as pool:
            for record, orfs in pool.imap(self._process_record, records):
                for j, orf in enumerate(orfs):
                    # wrap the protein into a Protein object
                    prot_seq = Seq(orf.translate())
                    protein = Protein(id=f"{record.id}_{j+1}", seq=prot_seq)
                    # wrap the gene into a Gene
                    yield Gene(
                        source=record,
                        start=min(orf.begin, orf.end),
                        end=max(orf.begin, orf.end),
                        strand=Strand.Coding if orf.strand == 1 else Strand.Reverse,
                        protein=protein,
                        qualifiers={
                            "inference": ["ab initio prediction:Prodigal:2.6"],
                            "transl_table": str(orf.translation_table),
                        }
                    )
