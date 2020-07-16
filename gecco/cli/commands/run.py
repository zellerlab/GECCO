"""Implementation of the ``gecco run`` subcommand.
"""

import glob
import logging
import multiprocessing
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
from ... import data
from ...data.hmms import Hmm, ForeignHmm
from ...hmmer import HMMER
from ...knn import ClusterKNN
from ...orf import PyrodigalFinder
from ...refine import ClusterRefiner


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
        -p <file>, --proteins <file>  a FASTA file containing proteins as
                                      input.

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
        --min-orfs <N>                how many ORFs are required for a
                                      sequence to be considered. [default: 5]
        -m <m>, --threshold <m>       the probability threshold for cluster
                                      detection. Default depends on the
                                      post-processing method (0.4 for gecco,
                                      0.6 for antismash).
        --postproc <method>           the method to use for cluster extraction
                                      (antismash or gecco). [default: gecco]

    Parameters - BGC Type Prediction:
        -d <d>, --distance <d>        the distance metric to use for kNN type
                                      prediction. [default: jensenshannon]
        -k <n>, --neighbors <n>       the number of neighbors to use for
                                      kNN type prediction [default: 5]

    Parameters - Debug:
        --model <model.crf>           the path to an alternative CRF model
                                      to use (obtained with `gecco train`).
        --hmm <lib.hmm.gz>            the path to one or more HMM libraries to
                                      use instead of the builtin ones.
    """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check value of numeric arguments
        self.args["--neighbors"] = int(self.args["--neighbors"])
        self.args["--min-orfs"] = int(self.args["--min-orfs"])
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
        self.args["--jobs"] = jobs = int(self.args["--jobs"])
        if jobs == 0:
            self.args["--jobs"] = multiprocessing.cpu_count()

        # Check the input exists
        input_ = self.args["--genome"] or self.args["--proteins"]
        if not os.path.exists(input_):
            self.logger.error("could not locate input file: {!r}", input)
            return 1

        # Check the hmms exist or use internal ones
        if self.args["--hmm"]:
            for hmm in self.args["--hmm"]:
                if not os.path.exists(hmm):
                    self.logger.error("could not locate hmm file: {!r}", hmm)
                    return 1
            self.args["--hmm"] = list(map(ForeignHmm, self.args["--hmm"]))
        else:
            self.args["--hmm"] = list(data.hmms.iter())

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
        format = guess_sequences_format(genome)
        sequences = list(SeqIO.parse(genome, format))

        self.logger.info("Findings genes in {} records", len(sequences))
        orf_finder = PyrodigalFinder(metagenome=True)
        genes = list(orf_finder.find_genes(sequences))
        self.logger.info("Found {} potential genes", len(genes))

        # --- HMMER ----------------------------------------------------------
        self.logger.info("Running domain annotation")

        # Run all HMMs over ORFs to annotate with protein domains
        def annotate(hmm: Hmm) -> "pandas.DataFrame":
            self.logger.debug("Starting annotation with HMM {} v{}", hmm.id, hmm.version)
            features = HMMER(hmm, self.args["--jobs"]).run(genes)
            self.logger.debug("Finished running HMM {}", hmm.id)

        with multiprocessing.pool.ThreadPool(self.args["--jobs"]) as pool:
            pool.map(annotate, self.args["--hmm"])

        # Count number of annotated domains
        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.logger.debug("Found {} domains across all proteins", count)

        # Filter i-evalue
        self.logger.debug("Filtering results with e-value under {}", self.args["--e-filter"])
        for gene in genes:
            gene.protein.domains = [
                d
                for d in gene.protein.domains
                if d.i_evalue < self.args["--e-filter"]
            ]
        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.logger.debug("Using remaining {} domains", count)

        # Sort genes
        self.logger.debug("Sort genes by coordinates")
        genes.sort(key=lambda g: (g.source.id, g.start, g.end))

        # --- CRF ------------------------------------------------------------
        self.logger.info("Predicting cluster probabilities with the CRF model")

        # Build the feature table from the list of annotated genes
        self.logger.debug("Building feature table from gene list")
        feats_df = pandas.concat([g.to_feature_table() for g in genes])
        feats_df.sort_values(by=["sequence_id", "start", "end", "domain_start"], inplace=True)

        # Load trained CRF model
        if self.args["--model"] is not None:
            self.logger.debug("Loading model from {!r}", self.args["--model"])
            with open(self.args["--model"], "rb") as bin:
                crf = pickle.load(bin)
        else:
            self.logger.debug("Loading model from package resources")
            crf = data.model.load()

        # Split input dataframe into one group per input sequence
        feats_df = crf.predict_marginals(data=[g for _, g in feats_df.groupby("sequence_id")])

        # Assign probabilities to data classes
        for gene in genes:
            rows = feats_df[(feats_df.protein_id == gene.id) & (feats_df.sequence_id == gene.source.id)]
            gene.probability = rows.p_pred.mean() if len(rows) else None

        # Write predictions to file
        pred_out = os.path.join(out_dir, f"{base}.features.tsv")
        self.logger.debug("Writing feature table to {!r}", pred_out)
        feats_df.to_csv(pred_out, sep="\t", index=False)

        # --- REFINE ---------------------------------------------------------
        self.logger.info("Extracting gene clusters from prediction")

        # Extract clusters from the predicted probability spectrum
        self.logger.debug("Using probability threshold of {}", self.args["--threshold"])
        refiner = ClusterRefiner(
             criterion=self.args["--postproc"],
            threshold=self.args["--threshold"],
            n_cds=self.args["--min-orfs"]
        )
        clusters = list(refiner.iter_clusters(genes))

        # Abort here if not clusters were found
        if clusters:
            self.logger.info("Found {} potential clusters", len(clusters))
        else:
            self.logger.warning("No gene clusters were found")
            return 0

        # # --- KNN ------------------------------------------------------------
        self.logger.info("Predicting BGC types")

        # Read embedded training matrix
        self.logger.debug("Reading embedded training matrix")
        training = data.knn.load_training_matrix()

        # Calculate new domain composition
        self.logger.debug("Calulating domain composition for each cluster")
        domains: Sequence[str] = training.domains # type: ignore
        new_comp = numpy.array([c.domain_composition(domains) for c in clusters])

        # Inititate kNN and predict types
        distance = self.args["--distance"]
        self.logger.debug("Running kNN classifier with {!r} metric", distance)
        knn = ClusterKNN(metric=distance, n_neighbors=self.args["--neighbors"])
        knn_pred = knn.fit_predict(training.compositions, new_comp, y=training.types)

        # Record predictions to the data classes
        for cluster, ty in zip(clusters, knn_pred):
            cluster.type, cluster.type_probability = ty

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
