import abc
import errno
import subprocess
import typing
from typing import Iterable, Type
from subprocess import DEVNULL

from ._meta import classproperty


class MissingBinaryError(RuntimeError):
    """A runtime error when a required program could not be found.
    """

    def __init__(self, name: str, args: Iterable[str]):
        self.name = name
        self.args = args


class BinaryRunner(metaclass=abc.ABCMeta):
    """A base class that wraps a binary program.
    """

    BINARY: typing.ClassVar[str] = NotImplemented

    @classmethod
    def has_binary(cls: Type["BinaryRunner"], args: Iterable[str] = ("--help",)) -> bool:
        try:
            cls.check_binary()
            return True
        except MissingBinaryError:
            return False

    @classmethod
    def check_binary(cls: Type["BinaryRunner"], args: Iterable[str] = ("--help",)) -> None:
        try:
            subprocess.run([cls.BINARY, *args], stdout=DEVNULL, stderr=DEVNULL)
        except OSError as err:
            if err.errno == errno.ENOENT:
                raise MissingBinaryError(cls.BINARY, args)
            raise

    def __init__(self):
        self.check_binary()
