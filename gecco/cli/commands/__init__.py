# PYTHON_ARGCOMPLETE_OK

import argparse
import functools
import signal
import typing
from typing import Optional, List, TextIO, Type, Callable, Iterable, Dict

from rich.console import Console

try:
    from rich_argparse import ArgumentDefaultsRichHelpFormatter
except ImportError:
    ArgumentDefaultsRichHelpFormatter = None

try:
    import argcomplete
except ImportError:
    argcomplete = None

from ... import __version__, __package__ as _program
from .._utils import patch_showwarnings
from .._log import showwarnings, make_logger, ConsoleLogger
from ._common import default_hmms
from ._parser import ConsoleHelpAction
from . import annotate, run, predict, train, cv, convert

if typing.TYPE_CHECKING:
    from ...hmmer import HMM
    from ...crf import ClusterCRF


def _formatter_class(console):
    if ArgumentDefaultsRichHelpFormatter is None:
        return argparse.ArgumentDefaultsHelpFormatter
    else:
        return functools.partial(ArgumentDefaultsRichHelpFormatter, console=console)


def configure_parser(
    parser: argparse.ArgumentParser,
    console: Console,
    program: str,
    version: str,
    *,
    defaults: Dict[str, object],
) -> argparse.ArgumentParser:
    _parser.configure_common(parser, console, program, version, main=True, defaults=defaults)

    commands = parser.add_subparsers(required=True, metavar="COMMAND")
    annotate.configure_parser(
        commands.add_parser(
            "annotate",
            formatter_class=_formatter_class(console),
            help="Annotate protein features of one or several contigs.",
            add_help=False,
        ),
        console,
        program,
        version,
        defaults=defaults,
    )
    run.configure_parser(
        commands.add_parser(
            "run",
            formatter_class=_formatter_class(console),
            help="Predict gene clusters from one or several contigs.",
            add_help=False,
        ),
        console,
        program,
        version,
        defaults=defaults,
    )
    predict.configure_parser(
        commands.add_parser(
            "predict",
            formatter_class=_formatter_class(console),
            help="Predict gene clusters on contigs that have been annotated.",
            add_help=False,
        ),
        console,
        program,
        version,
        defaults=defaults,
    )
    train.configure_parser(
        commands.add_parser(
            "train",
            formatter_class=_formatter_class(console),
            help="Train a new CRF model on pre-generated tables.",
            add_help=False,
        ),
        console,
        program,
        version,
        defaults=defaults,
    )
    cv.configure_parser(
        commands.add_parser(
            "cv",
            formatter_class=_formatter_class(console),
            help="Train and evaluate a model using cross-validation.",
            add_help=False,
        ),
        console,
        program,
        version,
        defaults=defaults,
    )
    convert.configure_parser(
        commands.add_parser(
            "convert",
            formatter_class=_formatter_class(console),
            help="Convert output files to a different format.",
            add_help=False,
        ),
        console,
        program,
        version,
        defaults=defaults,
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
    defaults: Optional[Dict[str, object]] = None,
) -> int:
    if defaults is None:
        defaults = {}

    parser = configure_parser(
        argparse.ArgumentParser(
            prog=program,
            formatter_class=_formatter_class(console),
            add_help=False,
        ),
        console,
        program,
        version,
        defaults=defaults,
    )
    if argcomplete is not None:
        argcomplete.autocomplete(parser)

    try:
        args = parser.parse_args(argv)
    except _parser.HelpExit:
        return 0

    if crf_type is None:
        from ...crf import ClusterCRF

        crf_type = ClusterCRF

    if classifier_type is None:
        from ...types import TypeClassifier

        classifier_type = TypeClassifier

    if console is None:
        console = Console(
            legacy_windows=not args.markup,
            no_color=not args.color,
            quiet=args.quiet,
            safe_box=not args.markup,
        )

    logger = make_logger(
        console,
        quiet=args.quiet,
        verbose=args.verbose,
        program=program,
    )

    with patch_showwarnings(
        functools.partial(showwarnings, logger, verbose=args.verbose, quiet=args.quiet)
    ):
        try:
            return args.run(
                args,
                logger,
                crf_type=crf_type,
                classifier_type=classifier_type,
                default_hmms=default_hmms,
            )
        except KeyboardInterrupt:
            logger.error("Interrupted")
            return -signal.SIGINT
        except OSError as e:
            logger.error(f"{e.strerror}: {e.filename!r}")
            return e.errno
        except RuntimeError as e:
            logger.error(f"{e.strerror}: {e.filename!r}")
            return 1
        except Exception as err:
            console.print_exception()
            if hasattr(err, "code"):
                return err.code
            elif hasattr(err, "errno"):
                return err.errno
            return 1
        else:
            return 0
