"""Simple data classes to expose embedded InterPro data.
"""

import gzip
import json
from dataclasses import dataclass, field, fields
from typing import Dict, List, Optional

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files  # type: ignore


__all__ = ["InterProEntry", "InterPro", "GeneOntologyTerm"]




@dataclass
class GOTerm:
    """A single term from the Gene Ontology.
    """
    accession: str
    name: str
    namespace: str

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
                self.accession == other.accession
            and self.name == other.name
            and self.namespace == other.namespace
        )

    def __hash__(self) -> int:
        return hash((type(self), self.accession, self.name, self.namespace))


@dataclass
class InterProEntry:
    """A single entry in the InterPro database.
    """

    accession: str
    members: List[str]
    name: str
    databases: List[str]
    type: str
    go_terms: List[GOTerm]
    go_functions: List[GOTerm]


@dataclass
class InterPro:
    """A subset of the InterPro database exposing domain metadata.
    """

    entries: List[InterProEntry]

    def __init__(self, entries: List[InterProEntry]):
        self.entries = entries
        self.by_accession = { member:entry for entry in entries for member in entry.members }

    @classmethod
    def load(cls) -> "InterPro":
        with files(__name__).joinpath("interpro.json").open() as f:
            data = json.load(f)
            entries = []
            for raw_entry in data:
                # get terms corresponding to domain
                go_terms = [
                    GOTerm(**t) for t in raw_entry.pop("go_terms")
                ]
                go_functions = [
                    GOTerm(**t, namespace="molecular_function") 
                    for t in raw_entry.pop("go_functions")
                ]
                entries.append(
                    InterProEntry(**raw_entry, go_terms=go_terms, go_functions=go_functions)
                )
        return cls(entries)
