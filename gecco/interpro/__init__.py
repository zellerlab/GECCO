"""Simple data classes to expose embedded InterPro data.
"""

import gzip
import json
from dataclasses import dataclass
from typing import Dict, List, Optional

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources  # type: ignore


__all__ = ["InterProEntry", "InterPro", "GeneOntologyTerm"]




@dataclass
class GeneOntologyTerm:
    """A single term from the Gene Ontology.
    """
    accession: str
    name: str
    namespace: str


@dataclass
class InterProEntry:
    """A single entry in the InterPro database.
    """

    accession: str
    members: List[str]
    name: str
    databases: List[str]
    type: str
    go_terms: List[GeneOntologyTerm]
    go_families: Dict[str, GeneOntologyTerm]


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
        with importlib_resources.open_binary(__name__, "interpro.json") as f:
            data = json.load(f)
            entries = []
            for raw_entry in data:
                go_terms = [
                    GeneOntologyTerm(**go_term)
                    for go_term in raw_entry.pop("go_terms")
                ]
                go_families = {
                    k:[
                        GeneOntologyTerm(namespace=k, **go_term)
                        for go_term in go_terms
                    ]
                    for k, go_terms in raw_entry.pop("go_families").items()
                }
                entries.append(
                    InterProEntry(**raw_entry, go_terms=go_terms, go_families=go_families)
                )
        return cls(entries)
