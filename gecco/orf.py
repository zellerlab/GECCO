"""Generic protocol for ORF detection in DNA sequences.
"""

import abc
import io
import os
import queue
import threading
import tempfile
import typing
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

    class _Worker(threading.Thread):

        @staticmethod
        def _none_callback(record, found, total):
            pass

        def __init__(
            self,
            metagenome: bool,
            record_count: "Value",
            record_queue: "queue.Queue[typing.Optional[SeqRecord]]",
            genes_queue: "queue.Queue[Gene]",
            callback: Optional[Callable[[SeqRecord, int], None]],
        ) -> None:
            super().__init__()
            self.pyrodigal = pyrodigal.Pyrodigal(meta=metagenome)
            self.record_queue = record_queue
            self.record_count = record_count
            self.genes_queue = genes_queue
            self.callback = callback or self._none_callback

        def run(self) -> None:
            while True:
                record = self.record_queue.get()
                if record is None:
                    self.record_queue.task_done()
                    return

                orfs = self.pyrodigal.find_genes(str(record.seq))
                for j, orf in enumerate(orfs):
                    # wrap the protein into a Protein object
                    prot_seq = Seq(orf.translate())
                    protein = Protein(id=f"{record.id}_{j+1}", seq=prot_seq)
                    # wrap the gene into a Gene
                    self.genes_queue.put(Gene(
                        source=record,
                        start=min(orf.begin, orf.end),
                        end=max(orf.begin, orf.end),
                        strand=Strand.Coding if orf.strand == 1 else Strand.Reverse,
                        protein=protein,
                        qualifiers={
                            "inference": ["ab initio prediction:Prodigal:2.6"],
                            "transl_table": str(orf.translation_table),
                        }
                    ))

                self.record_queue.task_done()
                self.callback(record, len(orfs), self.record_count.value)

    def __init__(self, metagenome: bool = True, cpus: int = 0) -> None:
        """Create a new `PyrodigalFinder` instance.

        Arguments:
            metagenome (bool): Whether or not to run PRODIGAL in metagenome
                mode, defaults to `True`.
            cpus (int): The number of threads to use to run Pyrodigal in
                parallel. Pass ``0`` to use the number of CPUs on the machine.

        """
        super().__init__()
        self.metagenome = metagenome
        self.cpus = cpus

    def find_genes(
        self,
        records: Iterable[SeqRecord],
        progress: Optional[Callable[[SeqRecord, int], None]] = None,
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

        """
        # count the number of CPUs to use
        _cpus = self.cpus if self.cpus > 0 else os.cpu_count()

        # create the queue to pass the objects around
        record_count = Value('i')
        record_queue = typing.cast("queue.Queue[typing.Optional[SeqRecord]]", queue.Queue())
        genes_queue = typing.cast("queue.Queue[SeqRecord]", queue.Queue())

        # create and launch one thread per CPU
        threads = []
        for _ in range(_cpus):
            thread = self._Worker(
                self.metagenome,
                record_count,
                record_queue,
                genes_queue,
                callback=progress,
            )
            thread.start()
            threads.append(thread)

        # queue the sequences passed as arguments
        for record in records:
            record_count.value += 1
            record_queue.put(record)

        # poison-pill the queue so that threads terminate when they
        # have consumed all the sequences
        for _ in threads:
            record_queue.put(None)

        # wait for all sequences to be processed
        record_queue.join()

        # return an iterable over the output
        sentinel = object()
        genes_queue.put(sentinel)
        return iter(genes_queue.get, sentinel)
