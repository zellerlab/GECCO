# PYTHON_ARGCOMPLETE_OK

import argparse
import functools
import signal
import typing
from typing import Optional, List, TextIO, Type, Callable, Iterable

from rich.console import Console

try:
    from rich_argparse import ArgumentDefaultsRichHelpFormatter as DefaultFormatter
except ImportError:
    from argparse import ArgumentDefaultsHelpFormatter as DefaultFormatter

try:
    import argcomplete
except ImportError:
    argcomplete = None

from ... import __version__, __package__ as _program
from .._utils import patch_showwarnings
from .._log import _showwarnings
from ._common import default_hmms
from . import annotate, run, predict, train, cv, convert

if typing.TYPE_CHECKING:
    from ...hmmer import HMM
    from ...crf import ClusterCRF


def configure_parser(
    parser: argparse.ArgumentParser,
    program: str,
    version: str,
) -> argparse.ArgumentParser:
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"{program} {version}",
        help="Show the program version number and exit.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        help="Increase the console output",
        default=0,
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="count",
        help="Reduce or disable the console output",
        default=0,
    )
    parser.add_argument(
        "--no-color",
        dest="color",
        action="store_false",
        help="Disable the console color",
    )
    parser.add_argument(
        "--no-markup",
        dest="markup",
        action="store_false",
        help="Disable the console markup",
    )

    commands = parser.add_subparsers(required=True, metavar="COMMAND")
    annotate.configure_parser(
        commands.add_parser(
            "annotate",
            formatter_class=DefaultFormatter,
            help="Annotate protein features of one or several contigs.",
            add_help=False,
        )
    )
    run.configure_parser(
        commands.add_parser(
            "run",
            formatter_class=DefaultFormatter,
            help="Predict gene clusters from one or several contigs.",
            add_help=False,
        )
    )
    predict.configure_parser(
        commands.add_parser(
            "predict",
            formatter_class=DefaultFormatter,
            help="Predict gene clusters on contigs that have been annotated.",
            add_help=False,
        )
    )
    train.configure_parser(
        commands.add_parser(
            "train",
            formatter_class=DefaultFormatter,
            help="Train a new CRF model on pre-generated tables.",
            add_help=False,
        )
    )
    cv.configure_parser(
        commands.add_parser(
            "cv",
            formatter_class=DefaultFormatter,
            help="Train and evaluate a model using cross-validation.",
            add_help=False,
        )
    )
    convert.configure_parser(
        commands.add_parser(
            "convert",
            formatter_class=DefaultFormatter,
            help="Convert output files to a different format.",
            add_help=False,
        )
    )

    return parser


def main(
    argv: Optional[List[str]] = None,
    console: Optional[Console] = None,
    *,
    program: str = _program,
    version: str = __version__,
    default_hmms: Callable[[], Iterable["HMM"]] = default_hmms,
    crf_type: Optional[Type["ClusterCRF"]] = None,
    classifier_type: Optional[Type["TypeClassifier"]] = None,
) -> int:
    if crf_type is None:
        from ...crf import ClusterCRF
        from ...types import TypeClassifier

        crf_type = ClusterCRF
        classifier_type = TypeClassifier

    parser = configure_parser(
        argparse.ArgumentParser(
            prog=program,
            formatter_class=DefaultFormatter,
            add_help=False,
        ),
        program,
        version,
    )
    if argcomplete is not None:
        argcomplete.autocomplete(parser)
    args = parser.parse_args(argv)

    if console is None:
        console = Console(
            legacy_windows=not args.markup,
            no_color=not args.color,
            quiet=args.quiet,
            safe_box=not args.markup,
        )

    with patch_showwarnings(
        functools.partial(
            _showwarnings, console, verbose=args.verbose, quiet=args.quiet
        )
    ):
        try:
            return args.run(
                args,
                console,
                crf_type=crf_type,
                classifier_type=classifier_type,
                default_hmms=default_hmms,
            )
        except Exception as err:
            console.print_exception()
            return getattr(err, "code", 1)
        else:
            return 0
