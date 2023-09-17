"""Generic protocol for ORF detection in DNA sequences.
"""

import abc
import io
import itertools
import os
import queue
import tempfile
import typing
from multiprocessing.pool import Pool, ThreadPool
from multiprocessing.sharedctypes import Value
from typing import Callable, Iterable, Iterator, List, Optional, Tuple, Type, Union

import Bio.SeqIO
import pyrodigal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .model import Gene, Protein, Strand


__all__ = ["ORFFinder", "PyrodigalFinder"]


class ORFFinder(metaclass=abc.ABCMeta):
    """An abstract base class to provide a generic ORF finder.
    """

    @abc.abstractmethod
    def find_genes(  # type: ignore
        self,
        records: Iterable[SeqRecord],
        progress: Optional[Callable[[SeqRecord, int], None]] = None,
    ) -> Iterable[Gene]:
        """Find all genes from a DNA sequence.
        """
        return NotImplemented


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
        self.orf_finder =  pyrodigal.GeneFinder(meta=metagenome, mask=mask)

    def _train(self, records: Iterable[SeqRecord]) -> pyrodigal.TrainingInfo:
        sequences = []
        for i, record in enumerate(records):
            if i > 0:
                sequences.append("TTAATTAATTAA")
            sequences.append(str(record.seq))
        if len(sequences) > 1:
            sequences.append("TTAATTAATTAA")
        return self.orf_finder.train("".join(sequences))

    def _process_record(self, record: SeqRecord) -> Tuple[SeqRecord, pyrodigal.Genes]:
        return record, self.orf_finder.find_genes(str(record.seq))

    def find_genes(
        self,
        records: Iterable[SeqRecord],
        progress: Optional[Callable[[SeqRecord, int], None]] = None,
        *,
        pool_factory: Union[Type[Pool], Callable[[Optional[int]], Pool]] = ThreadPool,
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
        _progress = (lambda x,y: None) if progress is None else progress

        # train first if needed
        if not self.metagenome:
            records = list(records)
            self._train(records)

        # run in parallel using a pool
        with pool_factory(_cpus) as pool:
            for record, orfs in pool.imap(self._process_record, records):
                _progress(record, len(orfs))
                for j, orf in enumerate(orfs):
                    # wrap the protein into a Protein object
                    prot_seq = Seq(orf.translate())
                    protein = Protein(id=f"{record.id}_{j+1}", seq=prot_seq)
                    # wrap the gene into a Gene
                    yield Gene(
                        source=record,
                        start=min(orf.begin, orf.end),
                        end=max(orf.begin, orf.end),
                        strand=Strand(orf.strand),
                        protein=protein,
                        qualifiers={
                            "inference": [f"ab initio prediction:Pyrodigal:{pyrodigal.__version__}"],
                            "transl_table": [str(orf.translation_table)],
                        }
                    )


class CDSFinder(ORFFinder):
    """An `ORFFinder` that simply extracts CDS annotations from records.
    """

    def __init__(
        self,
        feature: str = "CDS",
        translation_table: int = 11,
        locus_tag: str = "locus_tag",
    ):
        self.feature = feature
        self.translation_table = translation_table
        self.locus_tag = locus_tag

    def find_genes(
        self,
        records: Iterable[SeqRecord],
        progress: Optional[Callable[[SeqRecord, int], None]] = None,
    ) -> Iterator[Gene]:
        """Find all genes contained in a sequence of DNA records.
        """
        ids = set()
        _progress = (lambda x,y: None) if progress is None else progress

        for record in records:
            genes_found = 0
            features = filter(lambda feat: feat.type == self.feature, record.features)
            for i, feature in enumerate(features):
                # get the gene translation
                tt = feature.qualifiers.get("transl_table", [self.translation_table])[0]
                if "translation" in feature.qualifiers:
                    prot_seq = Seq(feature.qualifiers["translation"][0])
                else:
                    prot_seq = feature.location.extract(record.seq).translate(table=tt)
                # get the gene name
                if self.locus_tag in feature.qualifiers:
                    protein = Protein(id=feature.qualifiers[self.locus_tag][0], seq=prot_seq)
                else:
                    protein = Protein(id=f"{record.id}_{i+1}", seq=prot_seq)
                # check IDs are unique
                if protein.id in ids:
                    raise ValueError(f"Duplicate gene identifier found in {record.id!r}: {protein.id!r}")
                ids.add(protein.id)
                # wrap the gene into a Gene
                yield Gene(
                    source=record,
                    start=feature.location.start + 1,
                    end=feature.location.end,
                    strand=Strand(feature.location.strand),
                    protein=protein,
                )
                genes_found += 1
            _progress(record, genes_found)
