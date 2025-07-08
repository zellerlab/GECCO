"""A module containing some utilities only relevant to the CLI."""

import contextlib
import io
import os
import typing
import warnings
from types import ModuleType, TracebackType
from typing import Callable, Union, Iterator, TextIO, Type, Optional

import rich.progress
from rich.progress import Column
from rich.text import Text

from .._meta import classproperty, zopen

if typing.TYPE_CHECKING:
    ShowWarning = Callable[
        [Union[Warning, str], Type[Warning], str, int, Optional[TextIO], Optional[str]],
        None,
    ]


class MofNCompleteUnitColumn(rich.progress.MofNCompleteColumn):

    def render(self, task: "Task") -> Text:
        completed = int(task.completed)
        total = int(task.total) if task.total is not None else "?"
        total_width = len(str(total))
        return Text(
            f"{completed:{total_width}d}{self.separator}{total} {task.fields['unit']}",
            style="progress.download",
        )


@contextlib.contextmanager
def patch_showwarnings(new_showwarning: "ShowWarning") -> Iterator[None]:
    """Make a context patching `warnings.showwarning` with the given function."""
    old_showwarning: ShowWarning = warnings.showwarning
    try:
        warnings.showwarning = new_showwarning  # type: ignore
        yield
    finally:
        warnings.showwarning = old_showwarning  # type: ignore


def guess_sequences_format(path: Union[str, "os.PathLike[str]"]) -> Optional[str]:
    """Guess the format of a sequence file located in ``path``.

    Supports the following formats:
        * GenBank
        * FASTA
        * EMBL

    """
    head = ""
    with zopen(path) as file:
        for head in iter(lambda: file.read(256), b""):
            head = head.strip()
            if head:
                break
    if head.startswith(b">"):
        return "fasta"
    elif head.startswith(b"LOCUS"):
        return "genbank"
    elif head.startswith(b"ID "):
        return "embl"
    else:
        return None
