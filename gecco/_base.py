import abc
import errno
import io
import subprocess
import typing
from typing import Iterable, Optional, Type, TextIO
from subprocess import DEVNULL

from ._meta import classproperty


class Dumpable(metaclass=abc.ABCMeta):
    """A metaclass for objects that can be dumped to a text file.
    """

    @abc.abstractmethod
    def dump(self, fh: TextIO) -> None:
        raise NotImplementedError

    def dumps(self) -> str:
        s = io.StringIO()
        self.dump(s)
        return s.getvalue()
