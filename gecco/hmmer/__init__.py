"""Compatibility wrapper for HMMER binaries and output.
"""
import collections
import configparser
import contextlib
import csv
import errno
import glob
import os
import re
import subprocess
import tempfile
import typing
from typing import Callable, Dict, Optional, Iterable, Iterator, List, Mapping, Type, Sequence

import pkg_resources
from Bio import SeqIO

from .._meta import requires
from .._base import BinaryRunner
from ..model import Gene, Domain
from ..interpro import InterPro

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

    _T = typing.TypeVar("_T", bound="DomainRow")


class DomainRow(typing.NamedTuple):
    """A single row in a domain table created by ``hmmsearch``.

    See Also:
        The description of each field in page 48 of the `HMMER manual
        <http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf>`_.

    """

    target_name: str
    target_accession: Optional[str]
    target_length: int
    query_name: str
    query_accession: Optional[str]
    query_length: int
    evalue: float
    score: float
    bias: float
    domain_number: int
    domain_total: int
    c_evalue: float
    i_evalue: float
    domain_score: float
    domain_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    post: float
    description: Optional[str]

    @classmethod
    def from_line(cls: Type["_T"], row: str) -> "_T":
        """Extract a `DomainRow` from a single domain table line.
        """
        line = list(filter(None, row.split(" ")))
        values = {}
        for i, (field, ty) in enumerate(cls.__annotations__.items()):
            if isinstance(ty, type):
                values[field] = ty(line[i])
            elif field == "target_accession" or field == "query_accession":
                values[field] = None if line[i] == "-" else line[i]
            elif field == "description":
                values[field] = " ".join(line[i:]) if line[i:] else None
        return cls(**values)


class HMM(typing.NamedTuple):
    """A Hidden Markov Model library to use with `~gecco.hmmer.HMMER`.
    """

    id: str
    version: str
    url: str
    path: str
    relabel_with: Optional[str] = None

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


class HMMER(BinaryRunner):
    """A wrapper for HMMER that scans a HMM library against protein sequences.
    """

    BINARY = "hmmsearch"

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

    def run(self, genes: Iterable[Gene]) -> List[Gene]:
        """Run HMMER on proteins of ``genes`` and update them with domains.

        Arguments:
            genes (iterable of `~gecco.model.Gene`): An iterable that yield
                genes to annotate with ``self.hmm``.

        """
        # collect genes and build an index of genes by protein id
        gene_index = collections.OrderedDict([(gene.id, gene) for gene in genes])

        # create a temporary file to write the input and output to
        seqs_tmp = tempfile.NamedTemporaryFile(prefix="hmmer", suffix=".faa")
        doms_tmp = tempfile.NamedTemporaryFile(prefix="hmmer", suffix=".dom", mode="rt")

        # write protein sequences
        protein_records = (g.protein.to_seq_record() for g in gene_index.values())
        SeqIO.write(protein_records, seqs_tmp.name, "fasta")

        # Prepare the command line arguments
        cmd = ["hmmsearch", "--noali", "--domtblout", doms_tmp.name]
        if self.cpus is not None:
            cmd.extend(["--cpu", str(self.cpus)])
        cmd.extend([self.hmm.path, seqs_tmp.name])

        # Run HMMER
        subprocess.run(cmd, stdout=subprocess.DEVNULL).check_returncode()

        # Load InterPro metadata for the annotation
        interpro = InterPro.load()

        # Read the domain table
        lines = filter(lambda line: not line.startswith("#"), doms_tmp)
        rows = map(DomainRow.from_line, lines)

        # update protein domains
        for row in rows:
            # extract domain from the domain table row
            accession = self.hmm.relabel(row.query_accession or row.query_name)
            entry = interpro.by_accession.get(accession)
            # add additional qualifiers with available metadata
            qualifiers: Dict[str, List[str]] = {
                "inference": ["protein motif"],
                "note": ["e-value: {}".format(row.i_evalue)],
                "db_xref": ["{}:{}".format(self.hmm.id.upper(), accession)],
                "function": [] if entry is None else [entry.name]
            }
            if entry is not None and entry.integrated is not None:
                qualifiers["db_xref"].append("InterPro:{}".format(entry.integrated))
            # add the domain to the protein domains of the right gene
            domain = Domain(accession, row.env_from, row.env_to, self.hmm.id, row.i_evalue, None, qualifiers)
            gene_index[row.target_name].protein.domains.append(domain)

        # return the updated list of genes that was given in argument
        return list(gene_index.values())


class PyHMMER(object):

    def __init__(self, hmm: HMM, cpus: Optional[int] = None) -> None:
        self.hmm = hmm
        self.cpus = cpus

    @requires("pyhmmer")
    def run(self, genes: Iterable[Gene], callback: Optional[Callable[..., None]]) -> List[Gene]:
        # collect genes and build an index of genes by protein id
        gene_index = collections.OrderedDict([(gene.id, gene) for gene in genes])

        # create a temporary file to write the input and output to
        seqs_tmp = tempfile.NamedTemporaryFile(prefix="hmmer", suffix=".faa")
        #doms_tmp = tempfile.NamedTemporaryFile(prefix="hmmer", suffix=".dom", mode="rt")

        # write protein sequences
        protein_records = (g.protein.to_seq_record() for g in gene_index.values())
        SeqIO.write(protein_records, seqs_tmp.name, "fasta")

        # parse easel sequences
        esl_abc = pyhmmer.easel.Alphabet.amino()
        esl_sqs = [ seq.digitize(esl_abc) for seq in pyhmmer.easel.SequenceFile(seqs_tmp.name)]

        # Run HMMER subprocess.run(cmd, stdout=subprocess.DEVNULL).check_returncode()
        with pyhmmer.plan7.HMMFile(self.hmm.path) as hmm_file:
            hits = pyhmmer.hmmsearch(hmm_file, esl_sqs, cpus=self.cpus, callback=callback)

        # Load InterPro metadata for the annotation
        interpro = InterPro.load()

        # Transcribe HMMER hits to GECCO model
        for hit in hits:
            target_name = hit.name.decode('utf-8')
            for domain in hit.domains:
                raw_acc = domain.alignment.hmm_name
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
                domain = Domain(accession, domain.env_from, domain.env_to, self.hmm.id, domain.i_evalue, None, qualifiers)
                gene_index[target_name].protein.domains.append(domain)

        # return the updated list of genes that was given in argument
        return list(gene_index.values())



def embedded_hmms() -> Iterator[HMM]:
    """Iterate over the embedded HMMs that are shipped with GECCO.
    """
    for ini in glob.glob(pkg_resources.resource_filename(__name__, "*.ini")):
        cfg = configparser.ConfigParser()
        cfg.read(ini)
        yield HMM(path=ini.replace(".ini", ".hmm"), **dict(cfg.items("hmm")))
