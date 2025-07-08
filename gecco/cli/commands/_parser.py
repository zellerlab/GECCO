import argparse
import pathlib


class HelpExit(SystemExit):
    pass


class ConsoleHelpAction(argparse.Action):

    def __init__(
        self,
        option_strings,
        dest=argparse.SUPPRESS,
        default=argparse.SUPPRESS,
        help=None,
        deprecated=False,
        console=None,
    ):
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help,
            deprecated=deprecated
        )
        self.console = console

    def __call__(self, parser, namespace, values, option_string=None):
        if self.console is None:
            parser.print_help()
        else:
            self.console.print(parser.format_help())
        raise HelpExit(0)



def configure_group_gene_calling(
    parser: argparse.ArgumentParser,
) -> "argparse.ArgumentGroup":
    group = parser.add_argument_group(
        "Gene Calling",
    )
    group.add_argument(
        "-M",
        "--mask",
        action="store_true",
        help=(
            "Enable unknown region masking to prevent genes from stretching "
            "across ambiguous nucleotides."
        ),
    )
    group.add_argument(
        "--cds-feature",
        type=str,
        help=(
            "Extract genes from annotated records using a feature rather "
            "than running de-novo gene-calling."
        ),
    )
    group.add_argument(
        "--locus-tag",
        type=str,
        default="locus_tag",
        help=(
            "The name of the feature qualifier to use for naming extracted "
            "genes when using the ``--cds-feature`` flag."
        ),
    )
    return group


def configure_group_domain_annotation(
    parser: argparse.ArgumentParser,
) -> "argparse.ArgumentGroup":
    group = parser.add_argument_group("Domain Annotation")
    group.add_argument(
        "--hmm",
        action="append",
        dest="hmms",
        type=pathlib.Path,
        help=(
            "The path to one or more alternative HMM file to use " "(in HMMER format)."
        ),
    )
    group.add_argument(
        "-e",
        "--e-filter",
        type=float,
        default=None,
        help=(
            "The e-value cutoff for protein domains to be included. This is "
            "not stable across versions, so consider using a p-value filter "
            "instead."
        ),
    )
    group.add_argument(
        "-p",
        "--p-filter",
        type=float,
        default=1e-9,
        help=("The p-value cutoff for protein domains to be included."),
    )
    group.add_argument(
        "--bit-cutoffs",
        choices=("noise", "gathering", "trusted"),
        help=("Use HMM-specific bitscore cutoffs to filter domain annotations."),
    )
    # --disentangle creates issues wrt. HMM annotation because it may lead to
    # different domains being selected and applied to a domain depending on
    # whether we annotate with the full HMM or domains from a whitelist only
    # so we disable and hide it by default.
    group.add_argument(
        "--disentangle",
        action="store_true",
        help=argparse.SUPPRESS,
        # help=(
        #     "Disentangle overlapping domains in each gene by keeping only "
        #     "the domains with the lowest E-value over a given position."
        # ),
    )
    return group


def configure_group_table_output(
    parser: argparse.ArgumentParser,
) -> "argparse.ArgumentGroup":
    group = parser.add_argument_group(
        "Output",
    )
    group.add_argument(
        "-o",
        "--output-dir",
        help="The directory in which to write the output files.",
        type=pathlib.Path,
        default=".",
    )
    group.add_argument(
        "--force-tsv",
        action="store_true",
        help=(
            "Always write TSV output files even when they are empty "
            "(e.g. because no genes or no clusters were found)."
        ),
    )
    return group


def configure_group_cluster_detection(
    parser: argparse.ArgumentParser,
) -> "argparse.ArgumentGroup":
    group = parser.add_argument_group(
        "Cluster Detection",
    )
    group.add_argument(
        "--model",
        type=pathlib.Path,
        help=(
            "The path to an alternative CRF model to use (obtained with "
            "`gecco train`)."
        ),
    )
    group.add_argument(
        "--no-pad",
        dest="pad",
        action="store_false",
        help=(
            "Disable padding of gene sequences (used to predict gene clusters "
            "in contigs smaller than the CRF window length)."
        ),
    )
    group.add_argument(
        "-c",
        "--cds",
        type=int,
        default=3,
        help=(
            "The minimum number of coding sequences a valid cluster must "
            "contain to be retained."
        ),
    )
    group.add_argument(
        "-m",
        "--threshold",
        type=float,
        default=0.8,
        help=("The probability threshold for cluster detection."),
    )
    group.add_argument(
        "--postproc",
        choices={"antismash", "gecco"},
        default="gecco",
        help=argparse.SUPPRESS,
        # help=(
        #     "The method to use for cluster validation."
        # )
    )
    group.add_argument(
        "-E",
        "--edge-distance",
        default=0,
        type=int,
        help=(
            "The minimum number of annotated genes that must separate a "
            "cluster from the edge. Edge clusters will still be included "
            "if they are longer. A lower number will increase the number "
            "of false positives on small contigs."
        ),
    )
    return group


def configure_group_training_data(
    parser: argparse.ArgumentParser,
) -> "argparse.ArgumentGroup":
    group = parser.add_argument_group("Training Data")
    group.add_argument(
        "--no-shuffle",
        action="store_false",
        dest="shuffle",
        help=("Disable shuffling the data before fitting the model."),
    )
    group.add_argument(
        "--seed",
        type=int,
        default=42,
        help=(
            "The seed to initialize the random number generator used for "
            "shuffling operations."
        ),
    )
    return group


def configure_group_training_parameters(
    parser: argparse.ArgumentParser,
) -> "argparse.ArgumentGroup":
    group = parser.add_argument_group("Training Parameters")
    group.add_argument(
        "-W",
        "--window-size",
        type=int,
        default=5,
        help=("The length of the sliding window for CRF predictions."),
    )
    group.add_argument(
        "--window-step",
        type=int,
        default=1,
        help=("The step of the sliding window for CRF predictions."),
    )
    group.add_argument(
        "--c1",
        type=float,
        default=0.15,
        help=("The strength of the L1 regularization."),
    )
    group.add_argument(
        "--c2",
        type=float,
        default=0.15,
        help=("The strength of the L2 regularization."),
    )
    group.add_argument(
        "--feature-type",
        choices=("protein", "domain"),
        default="protein",
        help=(
            "The level at which the features should be extracted and "
            "given to the CRF."
        ),
    )
    group.add_argument(
        "--select",
        type=float,
        default=None,
        help=(
            "The fraction of most significant features to select from "
            "the training data prior to training the CRF."
        ),
    )
    group.add_argument(
        "--correction",
        type=str,
        default=None,
        help=(
            "The multiple test correction method to use when computing "
            "significance with multiple testing."
        ),
    )
    return group


def configure_group_input_tables(
    parser: argparse.ArgumentParser,
) -> "argparse.ArgumentGroup":
    group = parser.add_argument_group(
        "Input Tables",
    )
    group.add_argument(
        "-f",
        "--features",
        type=pathlib.Path,
        action="append",
        required=True,
        help=("The path to a domain annotation table, used to train the CRF model."),
    )
    group.add_argument(
        "-g",
        "--genes",
        type=pathlib.Path,
        required=True,
        help=(
            "The path to a gene table, containing the coordinates of the genes "
            "inside the training sequence."
        ),
    )
    group.add_argument(
        "-c",
        "--clusters",
        type=pathlib.Path,
        required=True,
        help=(
            "The path to a cluster annotation table, used to extract the "
            "domain composition for the type classifier."
        ),
    )
    return group
