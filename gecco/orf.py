"""Generic protocol for ORF detection in DNA sequences.
"""

import abc
import io
import itertools
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
            record_queue: "queue.Queue[typing.Optional[typing.Tuple[int, SeqRecord]]]",
            orfs_queue: "queue.Queue[Gene]",
            kill_switch: threading.Event,
            orfs_found: typing.List[threading.Event],
            callback: Optional[Callable[[SeqRecord, int], None]],
        ) -> None:
            super().__init__()
            self.pyrodigal = pyrodigal.Pyrodigal(meta=metagenome)
            self.record_queue = record_queue
            self.record_count = record_count
            self.orfs_queue = orfs_queue
            self.orfs_found = orfs_found
            self.callback = callback or self._none_callback
            self.error: typing.Optional[BaseException] = None
            self.kill_switch = kill_switch

        def kill(self) -> None:
            self.kill_switch.set()

        def run(self) -> None:
            while not self.kill_switch.is_set():
                args = self.record_queue.get()
                if args is None:
                    self.record_queue.task_done()
                    return
                else:
                    index, record = args
                try:
                    self.process(index, record)
                    self.record_queue.task_done()
                except BaseException as exc:
                    self.error = exc
                    self.kill()
                    return
                finally:
                    self.orfs_found[index].set()

        def process(self, index, record):
            orfs = self.pyrodigal.find_genes(str(record.seq))
            genes = []
            for j, orf in enumerate(orfs):
                # wrap the protein into a Protein object
                prot_seq = Seq(orf.translate())
                protein = Protein(id=f"{record.id}_{j+1}", seq=prot_seq)
                # wrap the gene into a Gene
                genes.append(Gene(
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
            self.orfs_queue.put((index, genes))
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

        Yields:
            `~gecco.model.Gene`: An iterator over all the genes found in
            the given records.

        """
        # count the number of CPUs to use
        _cpus = self.cpus if self.cpus > 0 else os.cpu_count()

        # create the queue to pass the objects around
        orfs_found: "typing.List[threading.Event]" = []
        record_count = Value('i')
        record_queue = typing.cast("queue.Queue[typing.Optional[typing.Tuple[int, SeqRecord]]]", queue.Queue())
        orfs_queue = typing.cast("queue.PriorityQueue[typing.Tuple[int, SeqRecord]]", queue.PriorityQueue())
        kill_switch = threading.Event()

        # create and launch one thread per CPU
        threads = []
        for _ in range(_cpus):
            thread = self._Worker(
                self.metagenome,
                record_count,
                record_queue,
                orfs_queue,
                kill_switch,
                orfs_found,
                callback=progress,
            )
            thread.start()
            threads.append(thread)

        # catch exceptions to kill threads in the background before exiting
        try:
            # enumerate queries, so that we know the index of each query
            # and we can yield the results in the same order
            queries = enumerate(records)
            # initially feed one query per thread so that they can start
            # working before we enter the main loop
            for index, record in itertools.islice(queries, _cpus):
                record_count.value += 1
                orfs_found.append(threading.Event())
                record_queue.put((index, record))
            # alternate between feeding queries to the threads and
            # yielding back results, if available
            orfs_yielded = 0
            while orfs_yielded < record_count.value:
                # get the next query record
                try:
                    index, query = next(queries)
                    record_count.value += 1
                    orfs_found.append(threading.Event())
                    record_queue.put((index, query))
                except StopIteration:
                    break
                # yield the top hits for the next query, if available
                if orfs_found[orfs_yielded].is_set():
                    yield from orfs_queue.get_nowait()[1]
                    orfs_yielded += 1
            # now that we exhausted all queries, poison pill the
            # threads so they stop on their own
            for _ in threads:
                record_queue.put(None)
            # yield remaining results
            while orfs_yielded < record_count.value:
                orfs_found[orfs_yielded].wait()
                yield from orfs_queue.get_nowait()[1]
                orfs_yielded += 1
        except queue.Empty:
            # the only way we can get queue.Empty is if a thread has set
            # the flag for `orfs_found[i]` without actually putting it in
            # the queue: this only happens when a background thread raised
            # an exception, so we must chain it
            for thread in threads:
                if thread.error is not None:
                    raise thread.error from None
        except BaseException:
            # make sure threads are killed to avoid being stuck,
            # e.g. after a KeyboardInterrupt
            kill_switch.set()
            raise
