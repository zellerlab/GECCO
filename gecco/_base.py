import abc
import csv
import errno
import io
import itertools
import operator
import subprocess
import typing
from collections.abc import Sized
from subprocess import DEVNULL
from typing import (
    Any,
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

from ._meta import classproperty, requires

if typing.TYPE_CHECKING:
    import pandas


_SELF = typing.TypeVar("_SELF")
_TABLE = typing.TypeVar("_TABLE", bound="Table")


def _parse_str(value: str) -> str:
    return value

def _parse_int(value: str) -> int:
    return int(value)

def _parse_float(value: str) -> float:
    return float(value)

def _parse_optional_float(value: str) -> typing.Optional[float]:
    return float(value) if value else None

def _parse_list_str(value: str) -> typing.List[str]:
    return value.split(";")

def _parse_optional_list_str(value: str) -> typing.Optional[typing.List[str]]:
    return None if not value else _parse_list_str(value)

def _format_int(value: int) -> str:
    return str(value)

def _format_str(value: str) -> str:
    return value

def _format_float(value: float) -> str:
    return str(value)

def _format_optional_float(value: typing.Optional[float]) -> str:
    return "" if value is None else str(value)

def _format_list_str(value: typing.List[str]) -> str:
    return ";".join(value)

def _format_optional_list_str(value: typing.Optional[typing.List[str]]) -> str:
    return "" if value is None else _format_list_str(value)


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


class Loadable(metaclass=abc.ABCMeta):
    """A metaclass for objects that can be loaded from a text file.
    """

    @classmethod
    @abc.abstractmethod
    def load(cls: typing.Type[_SELF], fh: TextIO) -> _SELF:
        raise NotImplementedError

    @classmethod
    def loads(cls: typing.Type[_SELF], s: str) -> _SELF:
        return cls.load(io.StringIO(s))  # type: ignore


class Table(Dumpable, Loadable, Sequence["Table.Row"]):
    """A metaclass for objects that
    """

    class Row(typing.NamedTuple):
        pass

    def __bool__(self) -> bool:  # noqa: D105
        return len(self) != 0

    def __iadd__(self: _TABLE, rhs: object) -> _TABLE:  # noqa: D105
        if not isinstance(rhs, type(self)):
            return NotImplemented
        for col in self.__annotations__:
            getattr(self, col).extend(getattr(rhs, col))
        return self

    @typing.overload
    def __getitem__(self, item: int) -> "Table.Row":  # noqa: D105
        pass

    @typing.overload
    def __getitem__(self: _TABLE, item: slice) -> _TABLE:  # noqa: D105
        pass

    def __getitem__(self: _TABLE, item: Union[slice, int]) -> Union[_TABLE, "Table.Row"]:   # noqa: D105
        columns = [getattr(self, col)[item] for col in self.__annotations__]
        if isinstance(item, slice):
            return type(self)(*columns)
        else:
            return self.Row(*columns)

    def __iter__(self) -> Iterator["Table.Row"]:  # noqa: D105
        columns = { c: operator.attrgetter(c) for c in self.__annotations__ }
        for i in range(len(self)):
            row = { c: getter(self)[i] for c, getter in columns.items() }
            yield self.Row(**row)

    @classmethod
    def _optional_columns(cls) -> typing.Set[str]:
        optional_columns = set()
        for name, ty in cls.Row.__annotations__.items():
            if ty == Optional[int] or ty == Optional[float] or ty == Optional[List[str]]:
                optional_columns.add(name)
        return optional_columns

    _FORMAT_FIELD: Dict[Any, Callable[[Any], str]] = {
        str: _format_str,
        int: _format_int,
        float: _format_float,
        typing.Optional[float]: _format_optional_float,
        typing.List[str]: _format_list_str,
        typing.Optional[typing.List[str]]: _format_optional_list_str,
    }

    def dump(self, fh: TextIO, dialect: str = "excel-tab", header: bool = True) -> None:
        """Write the table in CSV format to the given file.

        Arguments:
            fh (file-like `object`): A writable file-handle opened in text mode
                to write the feature table to.
            dialect (`str`): The CSV dialect to use. See `csv.list_dialects`
                for allowed values.
            header (`bool`): Whether or not to include the column header when
                writing the table (useful for appending to an existing table).
                Defaults to `True`.

        """
        writer = csv.writer(fh, dialect=dialect)
        column_names = list(self.__annotations__)
        optional = self._optional_columns()

        # do not write optional columns if they are completely empty
        for name in optional:
            if all(x is None for x in getattr(self, name)):
                column_names.remove(name)

        # write header if desired
        if header:
            writer.writerow(column_names)

        # write each row
        columns = [getattr(self, name) for name in column_names]
        formatters = [self._FORMAT_FIELD[self.Row.__annotations__[name]] for name in column_names]
        for i in range(len(self)):
            writer.writerow([format(col[i]) for col,format in zip(columns, formatters)])

    _PARSE_FIELD: Dict[Any, Callable[[str], Any]] = {
        str: _parse_str,
        int: _parse_int,
        float: _parse_float,
        typing.Optional[float]: _parse_optional_float,
        typing.List[str]: _parse_list_str,
        typing.Optional[typing.List[str]]: _parse_optional_list_str,
    }

    @classmethod
    def load(cls: typing.Type[_TABLE], fh: TextIO, dialect: str = "excel-tab") -> _TABLE:
        """Load a table in CSV format from a file handle in text mode.
        """
        table = cls()
        reader = csv.reader(fh, dialect=dialect)
        header = next(reader)

        # get the name of each column and check which columns are optional
        columns = [getattr(table, col) for col in header]
        optional = cls._optional_columns()
        parsers = [cls._PARSE_FIELD[table.Row.__annotations__[col]] for col in header]

        # check that if a column is missing, it is one of the optional values
        missing = set(cls.__annotations__).difference(header)
        missing_required = missing.difference(optional)
        if missing_required:
            raise ValueError("table is missing columns: {}".format(", ".join(missing_required)))

        # extract elements from the CSV rows
        for row in reader:
            for col, value, parse in itertools.zip_longest(columns, row, parsers):
                col.append(parse(value))
        for col in missing:
            getattr(table, col).extend(None for _ in range(len(table)))
        return table

    @requires("pandas")
    def to_dataframe(self) -> "pandas.DataFrame":  # type: ignore
        """Convert the table to a `~pandas.DataFrame`.

        Raises:
            ImportError: if the `pandas` module could not be imported.

        """
        frame = pandas.DataFrame()
        for column in self.__annotations__:
            frame[column] = getattr(self, column)
        return frame
