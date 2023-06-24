"""Compatibility wrapper for HMMER binaries and output.
"""
import abc
import atexit
import collections
import configparser
import contextlib
import csv
import errno
import glob
import gzip
import itertools
import os
import re
import subprocess
import tempfile
import typing
from typing import Any, BinaryIO, Callable, Container, Dict, Optional, Iterable, Iterator, List, Mapping, Type, Sequence

import pyhmmer
from Bio import SeqIO
from pyhmmer.hmmer import hmmsearch

from .._meta import requires, UniversalContainer, zopen
from ..model import Gene, Domain
from ..interpro import InterPro

try:
    from importlib.resources import files, as_file
except ImportError:
    from importlib_resources import files, as_file  # type: ignore

__all__ = ["DomainAnnotator", "HMM", "PyHMMER", "embedded_hmms"]


class HMM(typing.NamedTuple):
    """A Hidden Markov Model library to use with `~gecco.hmmer.HMMER`.
    """

    id: str
    version: str
    url: str
    path: str
    size: Optional[int] = None
    relabel_with: Optional[str] = None
    md5: Optional[str] = None

    def relabel(self, domain: str) -> str:
        """Rename a domain obtained by this HMM to the right accession.

        This method can be used with HMM libraries that have separate HMMs
        for the same domain, such as Pfam.
        """
        if self.relabel_with is None:
            return domain
        before, after = re.match("^s/(.*)/(.*)/$", self.relabel_with).groups()  # type: ignore
        regex = re.compile(before)
        return regex.sub(after, domain)


class DomainAnnotator(metaclass=abc.ABCMeta):
    """An abstract class for annotating genes with protein domains.
    """

    def __init__(self, hmm: HMM, cpus: Optional[int] = None, whitelist: Optional[Container[str]] = None) -> None:
        """Prepare a new HMMER annotation handler with the given ``hmms``.

        Arguments:
            hmm (str): The path to the file containing the HMMs.
            cpus (int, optional): The number of CPUs to allocate for the
                ``hmmsearch`` command. Give ``None`` to use the default.
            whitelist (container of str): If given, a container containing
                the accessions of the individual HMMs to annotate with. If
                `None` is given, annotate with the entire file.

        """
        super().__init__()
        self.hmm = hmm
        self.cpus = cpus
        self.whitelist = UniversalContainer() if whitelist is None else whitelist

    @abc.abstractmethod
    def run(self, genes: Iterable[Gene]) -> List[Gene]:
        """Run annotation on proteins of ``genes`` and update their domains.

        Arguments:
            genes (iterable of `~gecco.model.Gene`): An iterable that yield
                genes to annotate with ``self.hmm``.

        """
        return NotImplemented


class PyHMMER(DomainAnnotator):
    """A domain annotator that uses `pyhmmer`.
    """

    def run(
        self,
        genes: Iterable[Gene],
        progress: Optional[Callable[[pyhmmer.plan7.HMM, int], None]] = None,
        bit_cutoffs: Optional[str] = None,
    ) -> List[Gene]:
        # collect genes and keep them in original order
        gene_index = list(genes)

        # convert proteins to Easel sequences, namind them after
        # their location in the original input to ignore any duplicate
        # protein identifiers
        esl_abc = pyhmmer.easel.Alphabet.amino()
        esl_sqs = [
            pyhmmer.easel.TextSequence(
                name=str(i).encode(),
                sequence=str(gene.protein.seq)
            ).digitize(esl_abc)
            for i, gene in enumerate(gene_index)
        ]

        with contextlib.ExitStack() as ctx:
            # decompress the input if needed
            file: BinaryIO = ctx.enter_context(zopen(self.hmm.path))
            # Only retain the HMMs which are in the whitelist
            hmm_file = ctx.enter_context(pyhmmer.plan7.HMMFile(file))
            profiles = (
                hmm
                for hmm in hmm_file
                if hmm.accession is None
                or self.hmm.relabel(hmm.accession.decode()) in self.whitelist
            )
            # Run search pipeline using the filtered HMMs
            cpus = 0 if self.cpus is None else self.cpus
            hmms_hits = hmmsearch(
                profiles,
                esl_sqs,
                cpus=cpus,
                callback=progress, # type: ignore
                Z=self.hmm.size,  # type: ignore
                domZ=self.hmm.size, # type: ignore
                bit_cutoffs=bit_cutoffs,  # type: ignore
            )

            # Load InterPro metadata for the annotation
            interpro = InterPro.load()

            # Transcribe HMMER hits to GECCO model
            for hits in hmms_hits:
                for hit in hits.reported:
                    target_index = int(hit.name)
                    for domain in hit.domains.reported:
                        # extract name and get InterPro metadata about hit
                        raw_acc = domain.alignment.hmm_accession or domain.alignment.hmm_name
                        accession = self.hmm.relabel(raw_acc.decode('utf-8'))
                        entry = interpro.by_accession.get(accession)

                        # extract coordinates
                        start = domain.alignment.target_from
                        end = domain.alignment.target_to

                        # extract qualifiers and GO terms
                        qualifiers: Dict[str, List[str]] = {
                            "inference": ["protein motif"],
                            "db_xref": ["{}:{}".format(self.hmm.id.upper(), accession)],
                            "note": [
                                "e-value: {}".format(domain.i_evalue),
                                "p-value: {}".format(domain.pvalue),
                            ],
                        }
                        if entry is not None:
                            qualifiers["function"] = [entry.name]
                            qualifiers["db_xref"].append("InterPro:{}".format(entry.accession))
                            go_terms = entry.go_terms
                            go_functions = entry.go_functions
                        else:
                            go_terms = []
                            go_functions = []

                        # add the domain to the protein domains of the right gene
                        assert domain.env_from < domain.env_to
                        assert domain.i_evalue >= 0
                        assert domain.pvalue >= 0
                        gene_index[target_index].protein.domains.append(
                            Domain(
                                accession,
                                start,
                                end,
                                self.hmm.id,
                                domain.i_evalue,
                                domain.pvalue,
                                go_terms=go_terms,
                                go_functions=go_functions,
                                qualifiers=qualifiers,
                            )
                        )

        # return the updated list of genes that was given in argument
        return gene_index


def embedded_hmms() -> Iterator[HMM]:
    """Iterate over the embedded HMMs that are shipped with GECCO.
    """
    for filename in files(__name__).glob("*.ini"):

        ini_ctx = as_file(filename)
        ini_path = ini_ctx.__enter__()
        atexit.register(ini_ctx.__exit__, None, None, None)

        cfg = configparser.ConfigParser()
        cfg.read(ini_path)
        args: Dict[str, Any] = dict(cfg.items("hmm"))
        size = int(args.pop("size", 0))

        hmm_ctx = as_file(filename.with_suffix(".h3m"))
        hmm_path = hmm_ctx.__enter__()
        atexit.register(hmm_ctx.__exit__, None, None, None)

        yield HMM(path=os.fspath(hmm_path), size=size, **args)
