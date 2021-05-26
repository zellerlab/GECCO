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
import itertools
import os
import re
import subprocess
import tempfile
import typing
from typing import Callable, Dict, Optional, Iterable, Iterator, List, Mapping, Type, Sequence

import pyhmmer
from Bio import SeqIO
from pyhmmer.hmmer import hmmsearch

from .._meta import requires
from ..model import Gene, Domain
from ..interpro import InterPro

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

__all__ = ["DomainAnnotator", "HMM", "PyHMMER", "embedded_hmms"]


class HMM(typing.NamedTuple):
    """A Hidden Markov Model library to use with `~gecco.hmmer.HMMER`.
    """

    id: str
    version: str
    url: str
    path: str
    size: int
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

    def __init__(self, hmm: HMM, cpus: Optional[int] = None) -> None:
        """Prepare a new HMMER annotation handler with the given ``hmms``.

        Arguments:
            hmm (str): The path to the file containing the HMMs.
            cpus (int, optional): The number of CPUs to allocate for the
                ``hmmsearch`` command. Give ``None`` to use the default.

        """
        super().__init__()
        self.hmm = hmm
        self.cpus = cpus

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
        progress: Optional[Callable[[pyhmmer.plan7.HMM, int], None]] = None
    ) -> List[Gene]:
        # collect genes and build an index of genes by protein id
        gene_index = collections.OrderedDict([(gene.id, gene) for gene in genes])

        # convert proteins to Easel sequences
        esl_abc = pyhmmer.easel.Alphabet.amino()
        esl_sqs = [
            pyhmmer.easel.TextSequence(
                name=gene.protein.id.encode(),
                sequence=str(gene.protein.seq)
            ).digitize(esl_abc)
            for gene in genes
        ]

        # Run search pipeline using the HMM
        with pyhmmer.plan7.HMMFile(self.hmm.path) as hmm_file:
            cpus = 0 if self.cpus is None else self.cpus
            hmms_hits = hmmsearch(
                hmm_file,
                esl_sqs,
                cpus=cpus,
                callback=progress,
                Z=len(esl_sqs),
                domZ=len(esl_sqs)
            )

            # Load InterPro metadata for the annotation
            interpro = InterPro.load()

            # Transcribe HMMER hits to GECCO model
            for hit in itertools.chain.from_iterable(hmms_hits):
                target_name = hit.name.decode('utf-8')
                for domain in hit.domains:
                    raw_acc = domain.alignment.hmm_accession or domain.alignment.hmm_name
                    accession = self.hmm.relabel(raw_acc.decode('utf-8'))
                    entry = interpro.by_accession.get(accession)

                    # extract qualifiers
                    qualifiers: Dict[str, List[str]] = {
                        "inference": ["protein motif"],
                        "note": ["e-value: {}".format(domain.i_evalue)],
                        "db_xref": ["{}:{}".format(self.hmm.id.upper(), accession)],
                        "function": [] if entry is None else [entry.name]
                    }
                    if entry is not None and entry.integrated is not None:
                        qualifiers["db_xref"].append("InterPro:{}".format(entry.integrated))

                    # add the domain to the protein domains of the right gene
                    assert domain.env_from < domain.env_to
                    assert domain.i_evalue >= 0
                    d = Domain(accession, domain.env_from, domain.env_to, self.hmm.id, domain.i_evalue, None, qualifiers)
                    gene_index[target_name].protein.domains.append(d)

        # return the updated list of genes that was given in argument
        return list(gene_index.values())


def embedded_hmms() -> Iterator[HMM]:
    """Iterate over the embedded HMMs that are shipped with GECCO.
    """
    for filename in importlib_resources.contents(__name__):

        if not filename.endswith(".ini"):
            continue

        ini_ctx = importlib_resources.path(__name__, filename)
        ini_path = ini_ctx.__enter__()
        atexit.register(ini_ctx.__exit__, None, None, None)

        cfg = configparser.ConfigParser()
        cfg.read(ini_path)
        args = dict(cfg.items("hmm"))
        args["size"] = int(args["size"])

        hmm_ctx = importlib_resources.path(__name__, filename.replace(".ini", ".h3m"))
        hmm_path = hmm_ctx.__enter__()
        atexit.register(hmm_ctx.__exit__, None, None, None)

        yield HMM(path=os.fspath(hmm_path), **args)
