"""Simple data classes to expose embedded InterPro data.
"""

import gzip
import json
from dataclasses import dataclass
from typing import Dict, List, Optional

import pkg_resources


@dataclass
class InterProEntry:
    """A single domain entry in the InterPro database.
    """

    accession: str
    go_terms: Dict[str, Dict[str, str]]
    integrated: Optional[str]
    member_databases: Dict[str, Dict[str, str]]
    name: str
    source_database: str
    type: str


@dataclass
class InterPro:
    """A subset of the InterPro database exposing domain metadata.
    """

    entries: List[InterProEntry]

    def __init__(self, entries: List[InterProEntry]):
        self.entries = entries
        self.by_accession = { entry.accession:entry for entry in entries }

    @classmethod
    def load(cls) -> "InterPro":
        with pkg_resources.resource_stream(__name__, "interpro.json.gz") as f:
            data = json.load(gzip.open(f, mode="rt"))
            entries = [ InterProEntry(**entry["metadata"]) for entry in data ]
        return cls(entries)
