"""Implementation of the ``gecco run`` subcommand.
"""

import glob
import itertools
import logging
import multiprocessing
import operator
import os
import pickle
import tempfile
import typing
import signal
from typing import Union

import numpy
import pandas
from Bio import SeqIO

from ._base import Command
from .._utils import guess_sequences_format
from ...crf import ClusterCRF
from ...hmmer import HMMER, embedded_hmms
from ...knn import ClusterKNN
from ...orf import PyrodigalFinder
from ...refine import ClusterRefiner
from ...model import Hmm


class Run(Command):  # noqa: D101

    summary = "predict BGC from one or several contigs."
    doc = f"""
    gecco run - {summary}

    Usage:
        gecco run (-h | --help)
        gecco run --genome <file> [--hmm <hmm>]... [options]

    Arguments:
        -g <file>, --genome <file>    a FASTA or GenBank file containing a
                                      genome as input.

    Parameters:
        -o <out>, --output-dir <out>  the directory in which to write the
                                      output files. [default: .]
        -j <jobs>, --jobs <jobs>      the number of CPUs to use for
                                      multithreading. Use 0 to use all of the
                                      available CPUs. [default: 0]

    Parameters - Domain Annotation:
        -e <e>, --e-filter <e>        the e-value cutoff for protein domains
                                      to be included. [default: 1e-5]

    Parameters - Cluster Detection:
        -c, --cds <N>                 the minimum number of coding sequences a
                                      valid cluster must contain. [default: 3]
        -m <m>, --threshold <m>       the probability threshold for cluster
                                      detection. Default depends on the
                                      post-processing method (0.4 for gecco,
                                      0.6 for antismash).
        --postproc <method>           the method to use for cluster validation
                                      (antismash or gecco). [default: gecco]

    Parameters - BGC Type Prediction:
        -d <d>, --distance <d>        the distance metric to use for kNN type
                                      prediction. [default: jensenshannon]
        -k <n>, --neighbors <n>       the number of neighbors to use for
                                      kNN type prediction [default: 5]

    Parameters - Debug:
        --model <model.crf>           the path to an alternative CRF model
                                      to use (obtained with `gecco train`).
    """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check value of numeric arguments
        self.args["--neighbors"] = int(self.args["--neighbors"])
        self.args["--cds"] = int(self.args["--cds"])
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1

        # Use default threshold value dependeing on postprocessing method
        if self.args["--threshold"] is None:
            if self.args["--postproc"] == "gecco":
                self.args["--threshold"] = 0.4
            elif self.args["--postproc"] == "antismash":
                self.args["--threshold"] = 0.6
        else:
            self.args["--threshold"] = float(self.args["--threshold"])

        # Check the `--cpu`flag
        self.args["--jobs"] = int(self.args["--jobs"]) or multiprocessing.cpu_count()

        # Check the input exists
        if not os.path.exists(self.args["--genome"]):
            self.logger.error("could not locate input file {!r}", input)
            return 1

        return None

    def __call__(self) -> int:  # noqa: D102
        # Make output directory
        out_dir = self.args["--output-dir"]
        self.logger.debug("Using output folder {!r}", out_dir)
        os.makedirs(out_dir, exist_ok=True)

        # --- ORFs -----------------------------------------------------------
        genome = self.args["--genome"]
        base, _ = os.path.splitext(os.path.basename(genome))

        self.logger.info("Loading sequences from genome file {!r}", genome)
        sequences = SeqIO.parse(genome, guess_sequences_format(genome))

        self.logger.debug("Extracting genes from input sequences")
        orf_finder = PyrodigalFinder(metagenome=True)
        genes = list(
            itertools.chain.from_iterable(map(orf_finder.find_genes, sequences))
        )
        self.logger.info("Found {} potential genes", len(genes))

        # --- HMMER ----------------------------------------------------------
        self.logger.info("Running domain annotation")

        # Run all HMMs over ORFs to annotate with protein domains
        def annotate(hmm: Hmm) -> "pandas.DataFrame":
            self.logger.debug(
                "Starting annotation with HMM {} v{}", hmm.id, hmm.version
            )
            features = HMMER(hmm, self.args["--jobs"]).run(genes)
            self.logger.debug("Finished running HMM {}", hmm.id)

        with multiprocessing.pool.ThreadPool(self.args["--jobs"]) as pool:
            pool.map(annotate, embedded_hmms())

        # Count number of annotated domains
        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.logger.debug("Found {} domains across all proteins", count)

        # Filter i-evalue
        self.logger.debug(
            "Filtering results with e-value under {}", self.args["--e-filter"]
        )
        for gene in genes:
            key = lambda d: d.i_evalue < self.args["--e-filter"]
            gene.protein.domains = list(filter(key, gene.protein.domains))

        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.logger.debug("Using remaining {} domains", count)

        # Sort genes
        self.logger.debug("Sort genes by coordinates")
        genes.sort(key=lambda g: (g.source.id, g.start, g.end))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

        # --- CRF ------------------------------------------------------------
        self.logger.info("Predicting cluster probabilities with the CRF model")

        self.logger.debug("Loading trained CRF model")
        crf = ClusterCRF.trained(self.args["--model"])

        self.logger.debug("Predicting BGC probabilities")
        genes = crf.predict_probabilities(genes)

        self.logger.debug("Extracting feature table")
        feats_df = pandas.concat([g.to_feature_table() for g in genes], sort=False)

        pred_out = os.path.join(out_dir, f"{base}.features.tsv")
        self.logger.debug("Writing feature table to {!r}", pred_out)
        feats_df.to_csv(pred_out, sep="\t", index=False)

        # --- REFINE ---------------------------------------------------------
        self.logger.info("Extracting gene clusters from prediction")

        # Extract clusters from the predicted probability spectrum
        self.logger.debug("Using probability threshold of {}", self.args["--threshold"])
        refiner = ClusterRefiner(
            self.args["--threshold"], self.args["--postproc"], self.args["--cds"]
        )
        clusters = list(refiner.iter_clusters(genes))

        # Abort here if not clusters were found
        if clusters:
            self.logger.info("Found {} potential clusters", len(clusters))
        else:
            self.logger.warning("No gene clusters were found")
            return 0

        # --- KNN ------------------------------------------------------------
        self.logger.info("Predicting BGC types")
        knn = ClusterKNN.trained(self.args["--model"], metric=self.args["--distance"])
        clusters = knn.predict_types(clusters)

        # --- RESULTS --------------------------------------------------------
        self.logger.info("Writing final result files to folder {!r}", out_dir)

        # Write predicted cluster coordinates to file
        cluster_out = os.path.join(out_dir, f"{base}.clusters.tsv")
        self.logger.debug("Writing cluster coordinates to {!r}", cluster_out)
        table = pandas.concat([c.to_cluster_table() for c in clusters])
        table.to_csv(cluster_out, sep="\t", index=False)

        # Write predicted cluster sequences to file
        for cluster in clusters:
            gbk_out = os.path.join(out_dir, f"{cluster.id}.gbk")
            self.logger.debug("Writing cluster {} to {!r}", cluster.id, gbk_out)
            SeqIO.write(cluster.to_record(), gbk_out, "genbank")

        # Exit gracefully
        self.logger.info("Successfully found {} clusters!", len(clusters))
        return 0
