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

    def __init__(self, name: str, params: Iterable[str]) -> None:
        self.params = list(params)
        self.name = name
        super().__init__(name, self.params)

    def __str__(self) -> str:
        return f"could not locate binary: {self.name}"


class BinaryRunner(metaclass=abc.ABCMeta):
    """A base class that wraps a binary program.
    """

    BINARY: typing.ClassVar[str] = NotImplemented

    @classmethod
    def has_binary(
        cls: Type["BinaryRunner"], args: Iterable[str] = ("--help",)
    ) -> bool:
        try:
            cls.check_binary()
            return True
        except MissingBinaryError:
            return False

    @classmethod
    def check_binary(
        cls: Type["BinaryRunner"], args: Iterable[str] = ("--help",)
    ) -> None:
        try:
            subprocess.run([cls.BINARY, *args], stdout=DEVNULL, stderr=DEVNULL)
        except OSError as err:
            if err.errno == errno.ENOENT:
                raise MissingBinaryError(cls.BINARY, args)
            raise

    def __init__(self) -> None:
        self.check_binary()
