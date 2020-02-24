import csv
import logging
import multiprocessing
import os
import pickle
import random
import typing

import numpy
import pandas
from Bio import SeqIO

from ._base import Command
from ... import data
from ...crf import ClusterCRF
from ...hmmer import HMMER
from ...orf import ORFFinder
from ...refine import ClusterRefiner
from ...preprocessing import truncate


class Annotate(Command):

    summary = "use HMMs to annotate features for some proteins."
    doc = f"""
    gecco annotate - {summary}

    Usage:
        gecco annotate (-h | --help)
        gecco annotate --genome <file> [--hmm <hmm>]... [options]
        gecco annotate --proteins <file> [--hmm <hmm>]...  [options]

    Arguments:
        -g <file>, --genome <file>    a FASTA or GenBank file containing a
                                      genome as input.
        -p <file>, --proteins <file>  a FASTA file containing proteins as
                                      input.

    Parameters:
        -o <out>, --output-dir <out>  the directory in which to write the
                                      output files. [default: .]
        -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                      multithreading. Use 0 to use all of the
                                      available CPUs. [default: 0]

    Parameters - Domain Annotation:
        --hmm <hmm>                   the path to an HMM file to use for
                                      domain annotation. Defaults to the
                                      internal PFam HMM.
        -e <e>, --e-filter <e>        the e-value cutoff for domains to
                                      be included. [default: 1e-5]
    """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check value of numeric arguments
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1

        # Check the `--cpu`flag
        self.args["--jobs"] = jobs = int(self.args["--jobs"])
        if jobs == 0:
            self.args["--jobs"] = multiprocessing.cpu_count()

        # Check the input exists
        input = self.args["--genome"] or self.args["--proteins"]
        if not os.path.exists(input):
            self.logger.error("could not locate input file: {!r}", input)
            return 1

        # Check the HMM file(s) exist.
        self.args["--hmm"] = self.args["--hmm"] or [data.realpath("hmms/Pfam-A.hmm.gz")]
        for hmm in self.args["--hmm"]:
            if not os.path.exists(hmm):
                self.logger.error("could not locate hmm file: {!r}", hmm)
                return 1

        return None

    def __call__(self) -> int:
        # Make output directory
        out_dir = self.args["--output-dir"]
        self.logger.debug("Using output folder: {!r}", out_dir)
        os.makedirs(out_dir, exist_ok=True)

        # --- ORFs -----------------------------------------------------------
        if self.args["--genome"] is not None:
            genome = self.args["--genome"]
            base, _ = os.path.splitext(os.path.basename(genome))

            prodigal_out = os.path.join(out_dir, "prodigal")
            self.logger.debug("Using PRODIGAL output folder: {!r}", prodigal_out)
            os.makedirs(prodigal_out, exist_ok=True)

            self.logger.info("Predicting ORFs with PRODIGAL")
            prodigal = ORFFinder(genome, prodigal_out, method="prodigal")
            orf_file = prodigal.run()
            prodigal = True

        else:
            orf_file = self.args["--proteins"]
            base, _ = os.path.splitext(os.path.basename(orf_file))
            prodigal = False

        # --- HMMER ----------------------------------------------------------
        self.logger.info("Running PFam domain annotation")

        # Run PFAM HMM DB over ORFs to annotate with Pfam domains
        features = []
        for hmm in self.args["--hmm"]:
            hmmer_out = os.path.join(out_dir, "hmmer", os.path.basename(hmm))
            os.makedirs(hmmer_out, exist_ok=True)
            hmmer = HMMER(orf_file, hmmer_out, hmm, prodigal, self.args["--jobs"])
            features.append(hmmer.run())

        feats_df = pandas.concat(features, ignore_index=True)
        self.logger.debug("Found {} domains across all proteins", len(feats_df))

        # Filter i-evalue
        self.logger.debug("Filtering results with e-value under {}", self.args["--e-filter"])
        feats_df = feats_df[feats_df["i_Evalue"] < self.args["--e-filter"]]
        self.logger.debug("Using remaining {} domains", len(feats_df))

        # Reformat pfam IDs
        feats_df = feats_df.assign(
            domain=feats_df["domain"].str.replace(r"(PF\d+)\.\d+", lambda m: m.group(1))
        )

        # Write feature table to file
        feat_out = os.path.join(out_dir, f"{base}.features.tsv")
        self.logger.debug("Writing feature table to {!r}", feat_out)
        feats_df.to_csv(feat_out, sep="\t", index=False)
