import abc
import csv
import errno
import io
import itertools
import math
import operator
import os
import subprocess
import typing
from collections.abc import Sized
from subprocess import DEVNULL
from typing import (
    Any,
    BinaryIO,
    Callable,
    Dict,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Type,
    TextIO,
    Union,
    Sequence,
)

import polars

try:
    _POLARS_VERSION = tuple(map(int, polars.__version__.split('.')))
except ValueError:
    _POLARS_VERSION = (0, 0, 0)

if typing.TYPE_CHECKING:
    _SELF = typing.TypeVar("_SELF")
    _TABLE = typing.TypeVar("_TABLE", bound="Table")


class Dumpable(metaclass=abc.ABCMeta):
    """A metaclass for objects that can be dumped to a file.
    """

    @abc.abstractmethod
    def dump(self, fh: Union[BinaryIO, str, os.PathLike]) -> None:
        raise NotImplementedError

    def dumps(self) -> bytes:
        s = io.BytesIO()
        self.dump(s)
        return s.getvalue()


class Loadable(metaclass=abc.ABCMeta):
    """A metaclass for objects that can be loaded from a file.
    """

    @classmethod
    @abc.abstractmethod
    def load(cls: typing.Type["_SELF"], fh: Union[BinaryIO, str, os.PathLike]) -> "_SELF":
        raise NotImplementedError

    @classmethod
    def loads(cls: typing.Type["_SELF"], s: bytes) -> "_SELF":
        return cls.load(io.BytesIO(s))  # type: ignore


class Table(Dumpable, Loadable):
    """A metaclass for objects that
    """

    class Column(typing.NamedTuple):
        name: str
        dtype: type
        default: Optional[object] = None

    @classmethod
    @abc.abstractmethod
    def _get_columns(cls) -> List["Table.Column"]:
        return []
    
    data: polars.DataFrame

    def __init__(self, data: Optional[polars.DataFrame] = None) -> None:
        columns = self._get_columns()

        if data is not None:
            for column in columns:
                if column.name not in data.columns:
                    data = data.with_columns(polars.lit(column.default).alias(column.name))
            self.data = data
        else:
            self.data = polars.DataFrame([
                polars.Series(column.name, dtype=column.dtype)
                for column in columns
            ])

    def __bool__(self) -> bool: # noqa: D105
        return len(self) != 0

    def __len__(self) -> int:
        return len(self.data)
    
    def __getattr__(self, name: str) -> object:
        try:
            return self.data[name]
        except polars.exceptions.ColumnNotFoundError as err:
            raise AttributeError(name) from err

    def __iadd__(self: "_TABLE", rhs: object) -> "_TABLE":  # noqa: D105
        if not isinstance(rhs, type(self)):
            return NotImplemented
        self.data = polars.concat([self.data, rhs.data])
        return self

    @classmethod
    def load(
        cls: typing.Type["_TABLE"], 
        fh: Union[BinaryIO, str, os.PathLike], 
    ) -> "_TABLE":
        columns = cls._get_columns()
        dtypes = { column.name: column.dtype for column in columns }
        if _POLARS_VERSION < (0, 16, 14):
            data = polars.read_csv(fh, sep="\t", dtypes=dtypes)
        else:
            data = polars.read_csv(fh, separator="\t", dtypes=dtypes)
        for column_name in data.columns:
            if data[column_name].dtype in (polars.Float32, polars.Float64):
                data = data.with_columns(polars.col(column_name).fill_null(math.nan))
        return cls(data)

    def dump(self, fh: Union[BinaryIO, str, os.PathLike]) -> None:
        # remove columns that contain only default values
        columns = [ column for column in self._get_columns() ]
        for column in columns.copy():
            if column.default is None:
                continue
            if math.isnan(column.default):
                if self.data[column.name].is_nan().all():
                    columns.remove(column)
            elif self.data[column.name].eq(column.default).all():
                columns.remove(column)
        # write the table as a TSV file
        view = self.data[[column.name for column in columns]]
        for column_name in view.columns:
            if view[column_name].dtype in (polars.Float32, polars.Float64):
                view = view.with_columns(polars.col(column_name).fill_nan(None))
        if _POLARS_VERSION < (0, 16, 14):
            view.write_csv(fh, sep="\t")
        else:
            view.write_csv(fh, separator="\t")
    