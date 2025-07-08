import argparse
import csv
import glob
import os
import pathlib
import re
from typing import Type, Callable, Iterable

import rich.progress
from rich.console import Console

try:
    from rich_argparse import ArgumentDefaultsRichHelpFormatter as DefaultFormatter
except ImportError:
    from argparse import ArgumentDefaultsHelpFormatter as DefaultFormatter

from .._log import ConsoleLogger


def configure_parser(parser: argparse.ArgumentParser):
    parser.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit."
    )

    commands = parser.add_subparsers(required=True, metavar="INPUT")

    gbk_parser = commands.add_parser(
        "gbk",
        formatter_class=DefaultFormatter,
        help="Convert the GenBank records to a different format.",
    )
    gbk_parser.set_defaults(input="gbk")
    gbk_parser.add_argument(
        "-i",
        "--input-dir",
        type=pathlib.Path,
        required=True,
        help=("The path to the input directory containing files to convert."),
    )
    gbk_parser.add_argument(
        "-f",
        "--format",
        required=True,
        choices=("bigslice", "fna", "faa"),
        help=("The output format to write."),
    )

    clusters_parser = commands.add_parser(
        "clusters",
        formatter_class=DefaultFormatter,
        help="Convert the clusters table to a different format.",
    )
    clusters_parser.set_defaults(input="clusters")
    clusters_parser.add_argument(
        "-i",
        "--input-dir",
        type=pathlib.Path,
        required=True,
        help=("The path to the input directory containing files to convert."),
    )
    clusters_parser.add_argument(
        "-f",
        "--format",
        required=True,
        choices=("gff",),
        help=("The output format to write."),
    )

    parser.set_defaults(run=run)

    return parser


# --- GenBank to BiG-SLICE GenBank ---------------------------------------------


def _convert_gbk_bigslice(
    logger: ConsoleLogger,
    input_dir: pathlib.Path,
    progress: rich.progress.Progress,
) -> None:
    import Bio.SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from ...model import ClusterTable

    # collect `*.clusters.tsv` files
    cluster_files = glob.glob(os.path.join(input_dir, "*.clusters.tsv"))
    unit = "table" if len(cluster_files) == 1 else "tables"
    task = progress.add_task(
        "Loading clusters", total=len(cluster_files), unit=unit, precision=""
    )
    # load the original coordinates from the `*.clusters.tsv` files
    coordinates = {}
    types = {}
    for cluster_file in progress.track(cluster_files, task_id=task):
        table = ClusterTable.load(cluster_file)
        for i, cluster_id in enumerate(table.cluster_id):
            coordinates[cluster_id] = (table.start[i], table.end[i])
            types[cluster_id] = table.type[i] or "Unknown"

    # collect `*_clusters_{N}.gbk` files
    gbk_files = glob.glob(os.path.join(input_dir, "*_cluster_*.gbk"))
    unit = "file" if len(gbk_files) == 1 else "files"
    task = progress.add_task(
        "Converting", total=len(gbk_files), unit=unit, precision=""
    )
    done = 0
    # rewrite GenBank files
    for gbk_file in progress.track(gbk_files, task_id=task, total=len(gbk_files)):
        # load record and ensure it comes from GECCO
        record = Bio.SeqIO.read(gbk_file, "genbank")
        if "GECCO-Data" not in record.annotations.get("structured_comment", {}):
            logger.warn(f"GenBank file {gbk_file!r} was not obtained by GECCO")
            continue
        # mark as AntiSMASH v5
        record.annotations["structured_comment"]["antiSMASH-Data"] = {
            "Version": "5.X",
            "Orig. start": coordinates[record.id][0],
            "Orig. end": coordinates[record.id][1],
        }
        # using a /subregion feature will make BiG-SLiCE think the
        # cluster is coming from MIBiG, and the predicted type will be
        # handled correctly
        subregion_feature = SeqFeature(
            FeatureLocation(0, len(record)), type="subregion"
        )
        subregion_feature.qualifiers["contig_edge"] = ["False"]
        subregion_feature.qualifiers["aStool"] = ["mibig"]
        subregion_feature.qualifiers["label"] = [types[record.id]]
        record.features.append(subregion_feature)
        # rename the {id}_cluster_{N}.gbk file to {id}.region{N}.gbk
        gbk_name = os.path.basename(gbk_file)
        contig_id, cluster_n = re.search(r"^(.*)_cluster_(\d+).gbk", gbk_name).groups()  # type: ignore
        new_name = os.path.join(
            input_dir, "{}.region{:03}.gbk".format(contig_id, int(cluster_n))
        )
        logger.info(f"Rewriting {gbk_file!r} to {new_name!r}")
        Bio.SeqIO.write(record, new_name, "genbank")
        done += 1
    logger.success("Converted", f"{done} GenBank {unit} to BiG-SLiCE format", level=0)


# --- GenBank to nucleotide FASTA ----------------------------------------------


def _convert_gbk_fna(
    logger: ConsoleLogger,
    input_dir: pathlib.Path,
    progress: rich.progress.Progress,
) -> None:
    import Bio.SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    # collect `*_clusters_{N}.gbk` files
    gbk_files = glob.glob(os.path.join(input_dir, "*_cluster_*.gbk"))
    unit = "file" if len(gbk_files) == 1 else "files"
    task = progress.add_task(
        "Converting", total=len(gbk_files), unit=unit, precision=""
    )
    done = 0
    # rewrite GenBank files
    for gbk_file in progress.track(gbk_files, task_id=task, total=len(gbk_files)):
        # load record and ensure it comes from GECCO
        record = Bio.SeqIO.read(gbk_file, "genbank")
        if "GECCO-Data" not in record.annotations.get("structured_comment", {}):
            logger.warn(f"GenBank file {gbk_file!r} was not obtained by GECCO")
            continue
        # convert to nucleotide FASTA
        new_name = re.sub(r"\.gbk$", ".fna", gbk_file)
        logger.info(f"Converting {gbk_file!r} to FASTA file {new_name!r}")
        Bio.SeqIO.write(record, new_name, "fasta")
        done += 1
    logger.success(
        "Converted", f"{done} GenBank {unit} to nucleotide FASTA format", level=0
    )


# --- GenBank to protein FASTA -------------------------------------------------


def _convert_gbk_faa(
    logger: ConsoleLogger,
    input_dir: pathlib.Path,
    progress: rich.progress.Progress,
) -> None:
    import Bio.SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    # collect `*_clusters_{N}.gbk` files
    gbk_files = glob.glob(os.path.join(input_dir, "*_cluster_*.gbk"))
    unit = "file" if len(gbk_files) == 1 else "files"
    task = progress.add_task(
        "Converting", total=len(gbk_files), unit=unit, precision=""
    )
    done = 0
    # rewrite GenBank files
    for gbk_file in progress.track(gbk_files, task_id=task, total=len(gbk_files)):
        # load record and ensure it comes from GECCO
        record = Bio.SeqIO.read(gbk_file, "genbank")
        if "GECCO-Data" not in record.annotations.get("structured_comment", {}):
            logger.warn(f"GenBank file {gbk_file!r} was not obtained by GECCO")
            continue
        # extract proteins
        proteins = []
        for feat in filter(lambda f: f.type == "CDS", record.features):
            if "locus_tag" not in feat.qualifiers:
                continue
            seq = Seq(feat.qualifiers["translation"][0])
            prot = SeqRecord(id=feat.qualifiers["locus_tag"][0], seq=seq)
            proteins.append(prot)
        # write proteins to FASTA
        new_name = re.sub(r"\.gbk$", ".faa", gbk_file)
        logger.info(f"Converting {gbk_file!r} proteins to FASTA file {new_name!r}")
        Bio.SeqIO.write(proteins, new_name, "fasta")
        done += 1
    logger.success(
        "Converted", f"{done} GenBank {unit} to protein FASTA format", level=0
    )


# --- Clusters TSV to GFF ------------------------------------------------------


def _convert_clusters_gff(
    logger: ConsoleLogger,
    input_dir: pathlib.Path,
    progress: rich.progress.Progress,
) -> None:
    import Bio.SeqIO
    from ...model import ClusterTable

    # collect `*_clusters_{N}.gbk` files
    tsv_files = glob.glob(os.path.join(input_dir, "*.clusters.tsv"))
    unit = "file" if len(tsv_files) == 1 else "files"
    task = progress.add_task(
        "Converting", total=len(tsv_files), unit=unit, precision=""
    )
    done = 0

    # rewrite GenBank files
    for tsv_file in progress.track(tsv_files, task_id=task, total=len(tsv_files)):
        table = ClusterTable.load(tsv_file)
        gff_file = re.sub(r"\.tsv$", ".gff", tsv_file)
        with open(gff_file, "w") as dst:
            writer = csv.writer(dst, dialect="excel-tab")
            writer.writerow(["##gff-version 3"])
            for row in table.data.rows(named=True):
                # load the GenBank files of the cluster to extract the GECCO version
                gbk_file = os.path.join(input_dir, f"{row['cluster_id']}.gbk")
                cluster = Bio.SeqIO.read(gbk_file, "genbank")
                annotations = cluster.annotations["structured_comment"]["GECCO-Data"]
                # make sure to have a BGC type to write down
                bgc_types = ["Unknown"] if not row["type"] else row["type"].split(";")
                # extract type probabilities
                type_probas = []
                for key, value in row.items():
                    if key.endswith("_probability"):
                        ty = key.split("_")[0].capitalize()
                        if ty == "Nrp":
                            ty = "NRP"
                        type_probas.append(f"Type{ty}={value}")
                # write the GFF row
                writer.writerow(
                    [
                        row["sequence_id"],
                        annotations["version"],
                        "BGC",
                        str(row["start"]),
                        str(row["end"]),
                        str(row["average_p"]),
                        ".",
                        ".",
                        ";".join(
                            [
                                f"ID={row['cluster_id']}",
                                f"Name={ '/'.join(sorted(bgc_types)) } cluster",
                                f"Type={ ','.join(sorted(bgc_types)) }",
                                f"ProbabilityAverage={row['average_p']}",
                                f"ProbabilityMax={row['max_p']}",
                                *type_probas,
                                f"Genes={row['proteins'].count(';') + 1}",
                                f"Domains={row['domains'].count(';') + 1}",
                            ]
                        ),
                    ]
                )
        done += 1
    logger.success("Converted", f"{done} TSV {unit} to GFF format", level=0)


# --- Run command --------------------------------------------------------------


def run(
    args: argparse.Namespace,
    console: Console,
    crf_type: Type["ClusterCRF"],
    default_hmms: Callable[[], Iterable["HMM"]],
) -> int:
    logger = ConsoleLogger(console, quiet=args.quiet, verbose=args.verbose)

    with rich.progress.Progress(console=logger.console) as progress:
        if args.input == "gbk":
            if args.format == "bigslice":
                _convert_gbk_bigslice(logger, args.input_dir, progress)
            elif args.format == "fna":
                _convert_gbk_fna(logger, args.input_dir, progress)
            elif args.format == "faa":
                _convert_gbk_faa(logger, args.input_dir, progress)
            else:
                raise NotImplementedError
        elif args.input == "clusters":
            if args.format == "gff":
                _convert_clusters_gff(logger, args.input_dir, progress)
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

    # exit successfully
    return 0
