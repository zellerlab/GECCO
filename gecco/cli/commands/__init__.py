# PYTHON_ARGCOMPLETE_OK

import argparse
import functools
import signal
import typing
from typing import Optional, List, TextIO, Type, Callable, Iterable

from rich.console import Console

# try:
from rich_argparse import ArgumentDefaultsRichHelpFormatter  # as DefaultFormatter

# except ImportError:
# from argparse import ArgumentDefaultsHelpFormatter as DefaultFormatter

try:
    import argcomplete
except ImportError:
    argcomplete = None

from ... import __version__, __package__ as _program
from .._utils import patch_showwarnings
from .._log import _showwarnings, ConsoleLogger
from ._common import default_hmms
from ._parser import ConsoleHelpAction
from . import annotate, run, predict, train, cv, convert

if typing.TYPE_CHECKING:
    from ...hmmer import HMM
    from ...crf import ClusterCRF


def configure_parser(
    parser: argparse.ArgumentParser,
    console: Console,
    program: str,
    version: str,
) -> argparse.ArgumentParser:
    parser.add_argument(
        "-h",
        "--help",
        action=ConsoleHelpAction,
        help="Show this help message and exit.",
        console=console,
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
            formatter_class=lambda prog: ArgumentDefaultsRichHelpFormatter(
                prog, console=console
            ),
            help="Annotate protein features of one or several contigs.",
            add_help=False,
        ),
        console,
    )
    run.configure_parser(
        commands.add_parser(
            "run",
            formatter_class=lambda prog: ArgumentDefaultsRichHelpFormatter(
                prog, console=console
            ),
            help="Predict gene clusters from one or several contigs.",
            add_help=False,
        ),
        console,
    )
    predict.configure_parser(
        commands.add_parser(
            "predict",
            formatter_class=lambda prog: ArgumentDefaultsRichHelpFormatter(
                prog, console=console
            ),
            help="Predict gene clusters on contigs that have been annotated.",
            add_help=False,
        ),
        console,
    )
    train.configure_parser(
        commands.add_parser(
            "train",
            formatter_class=lambda prog: ArgumentDefaultsRichHelpFormatter(
                prog, console=console
            ),
            help="Train a new CRF model on pre-generated tables.",
            add_help=False,
        ),
        console,
    )
    cv.configure_parser(
        commands.add_parser(
            "cv",
            formatter_class=lambda prog: ArgumentDefaultsRichHelpFormatter(
                prog, console=console
            ),
            help="Train and evaluate a model using cross-validation.",
            add_help=False,
        ),
        console,
    )
    convert.configure_parser(
        commands.add_parser(
            "convert",
            formatter_class=lambda prog: ArgumentDefaultsRichHelpFormatter(
                prog, console=console
            ),
            help="Convert output files to a different format.",
            add_help=False,
        ),
        console,
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

        crf_type = ClusterCRF

    if classifier_type is None:
        from ...types import TypeClassifier

        classifier_type = TypeClassifier

    parser = configure_parser(
        argparse.ArgumentParser(
            prog=program,
            formatter_class=lambda prog: ArgumentDefaultsRichHelpFormatter(
                prog, console=console
            ),
            add_help=False,
        ),
        console,
        program,
        version,
    )
    if argcomplete is not None:
        argcomplete.autocomplete(parser)

    try:
        args = parser.parse_args(argv)
    except _parser.HelpExit:
        return 0

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
        logger = ConsoleLogger(
            console,
            quiet=args.quiet,
            verbose=args.verbose,
            program=program,
        )
        try:
            return args.run(
                args,
                logger,
                crf_type=crf_type,
                classifier_type=classifier_type,
                default_hmms=default_hmms,
            )
        except Exception as err:
            console.print_exception()
            return getattr(err, "code", 1)
        else:
            return 0
