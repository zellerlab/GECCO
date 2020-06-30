"""Implementation of the ``gecco run`` subcommand.
"""

import csv
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

    summary = "predict BGC from a genome or from individual proteins."
    doc = f"""
    gecco run - {summary}

    Usage:
        gecco run (-h | --help)
        gecco run --genome <file>   [--hmm <hmm>]... [options]
        gecco run --proteins <file> [--hmm <hmm>]... [options]

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
        -e <e>, --e-filter <e>        the e-value cutoff for PFam domains to
                                      be included. [default: 1e-5]

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
        self.logger.debug("Using output folder: {!r}", out_dir)
        os.makedirs(out_dir, exist_ok=True)

        # --- ORFs -----------------------------------------------------------
        if self.args["--genome"] is not None:
            genome = self.args["--genome"]
            base, _ = os.path.splitext(os.path.basename(genome))

            self.logger.info("Loading sequences from genome: {!r}", genome)
            format = guess_sequences_format(genome)
            sequences = SeqIO.parse(genome, format)
            self.logger.info("Predicting ORFs with PRODIGAL")
            orf_finder = PyrodigalFinder(metagenome=True)
            proteins = orf_finder.find_proteins(sequences)

            # we need to keep all the ORFs in a file because we will need
            # them when extracting cluster sequences
            _orf_temp = tempfile.NamedTemporaryFile(prefix="gecco", suffix=".faa")
            orf_file = _orf_temp.name
            SeqIO.write(proteins, orf_file, "fasta")
            prodigal = True

        else:
            _orf_temp = None
            orf_file = self.args["--proteins"]
            base, _ = os.path.splitext(os.path.basename(orf_file))
            prodigal = False

        # count the number of detected proteins without keeping them all
        # in memory
        orf_index = SeqIO.index(orf_file, "fasta")
        self.logger.info("Found {} potential proteins", len(orf_index))

        # --- HMMER ----------------------------------------------------------
        self.logger.info("Running domain annotation")

        # Run all HMMs over ORFs to annotate with protein domains
        def annotate(hmm: Union[Hmm, ForeignHmm]) -> "pandas.DataFrame":
            self.logger.debug(
                "Starting annotation with HMM {} v{}", hmm.id, hmm.version
            )
            hmmer_out = os.path.join(out_dir, "hmmer", hmm.id)
            os.makedirs(hmmer_out, exist_ok=True)
            hmmer = HMMER(hmm.path, self.args["--jobs"])
            result = hmmer.run(SeqIO.parse(orf_file, "fasta"), prodigal=prodigal)
            self.logger.debug("Finished running HMM {}", hmm.id)
            return result.assign(hmm=hmm.id, domain=hmm.relabel(result.domain))

        with multiprocessing.pool.ThreadPool(self.args["--jobs"]) as pool:
            features = pool.map(annotate, self.args["--hmm"])

        feats_df = pandas.concat(features, ignore_index=True)
        self.logger.debug("Found {} domains across all proteins", len(feats_df))

        # Filter i-evalue
        self.logger.debug(
            "Filtering results with e-value under {}", self.args["--e-filter"]
        )
        feats_df = feats_df[feats_df["i_Evalue"] < self.args["--e-filter"]]
        self.logger.debug("Using remaining {} domains", len(feats_df))

        # Write feature table to file
        feat_out = os.path.join(out_dir, f"{base}.features.tsv")
        self.logger.debug("Writing feature table to {!r}", feat_out)
        feats_df.to_csv(feat_out, sep="\t", index=False)

        # --- CRF ------------------------------------------------------------
        self.logger.info("Predicting cluster probabilities with the CRF model")
        if self.args["--model"] is not None:
            with open(self.args["--model"], "rb") as bin:
                crf = pickle.load(bin)
        else:
            crf = data.load("model/crf.model")

        # Compute reverse i_Evalue to be used as weight
        feats_df["rev_i_Evalue"] = 1 - feats_df.i_Evalue

        # If extracted from genome, split input dataframe into sequence
        feats_df = crf.predict_marginals(
            data=[seq for _, seq in feats_df.groupby("sequence_id")]
        )

        # Write predictions to file
        pred_out = os.path.join(out_dir, f"{base}.pred.tsv")
        self.logger.debug("Writing cluster probabilities to {!r}", pred_out)
        feats_df.to_csv(pred_out, sep="\t", index=False)

        # --- REFINE ---------------------------------------------------------
        self.logger.info("Extracting clusters")
        self.logger.debug("Using probability threshold of {}", self.args["--threshold"])
        refiner = ClusterRefiner(threshold=self.args["--threshold"])

        clusters = []
        for sid, subdf in feats_df.groupby("sequence_id"):
            if len(subdf["protein_id"].unique()) < self.args["--min-orfs"]:
                self.logger.warn("Skipping sequence {!r} because it is too short", sid)
                continue
            clusters.extend(
                refiner.iter_clusters(
                    subdf, criterion=self.args["--postproc"], prefix=sid,
                )
            )

        if not clusters:
            self.logger.warning("No gene clusters were found")
            return 0

        # --- KNN ------------------------------------------------------------
        self.logger.info("Predicting BGC types")

        # Read embedded training matrix
        self.logger.debug("Reading embedded training matrix")
        training_matrix = data.realpath("knn/domain_composition.tsv.gz")
        train_df = pandas.read_csv(training_matrix, sep="\t", encoding="utf-8")
        train_comp = train_df.iloc[:, 2:].values
        id_array = train_df["BGC_id"].values
        types_array = train_df["BGC_type"]
        domains_array = train_df.columns.values[2:]

        # Calculate new domain composition
        self.logger.debug("Calulating domain composition for each cluster")
        new_comp = numpy.array([c.domain_composition(domains_array) for c in clusters])

        # Inititate kNN and predict types
        distance = self.args["--distance"]
        self.logger.debug("Running kNN classifier with metric: {!r}", distance)
        knn = ClusterKNN(metric=distance, n_neighbors=self.args["--neighbors"])
        knn_pred = knn.fit_predict(train_comp, new_comp, y=types_array)

        # --- RESULTS --------------------------------------------------------
        self.logger.info("Writing final results file")

        # Write predicted cluster coordinates to file
        cluster_out = os.path.join(out_dir, f"{base}.clusters.tsv")
        self.logger.debug("Writing cluster coordinates to {!r}", cluster_out)
        with open(cluster_out, "wt") as f:
            csv.writer(f, dialect="excel-tab").writerow(
                [
                    "sequence_id",
                    "BGC_id",
                    "start",
                    "end",
                    "average_p",
                    "max_p",
                    "BGC_type",
                    "BGC_type_p",
                    "proteins",
                    "domains",
                ]
            )
            for cluster, ty in zip(clusters, knn_pred):
                cluster.bgc_type, cluster.type_prob = ty
                cluster.write_to_file(f, long=True)

        # Write predicted cluster sequences to file
        for cluster in clusters:
            prots_out = os.path.join(out_dir, f"{cluster.name}.proteins.faa")
            self.logger.debug("Writing proteins of {} to {!r}", cluster.name, prots_out)
            with open(prots_out, "w") as out:
                for id_ in cluster.prot_ids:
                    p = orf_index[id_]
                    p.description = f"{cluster.name} # {p.description}"
                    SeqIO.write(p, out, "fasta")

        # Remove the temporary protein file is needed to get rid of resource
        # warnings
        if _orf_temp is not None:
            _orf_temp.close()

        # Exit gracefully
        self.logger.info("Successfully found {} clusters!", len(clusters))
        return 0
