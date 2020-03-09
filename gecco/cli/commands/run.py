import csv
import glob
import logging
import multiprocessing
import os
import pickle
import typing

import numpy
import pandas
from Bio import SeqIO

from ._base import Command
from ... import data
from ...crf import ClusterCRF
from ...hmmer import HMMER
from ...knn import ClusterKNN
from ...orf import ORFFinder
from ...refine import ClusterRefiner


class Run(Command):

    summary = "predict BGC from a genome or from individual proteins."
    doc = f"""
    gecco run - {summary}

    Usage:
        gecco run (-h | --help)
        gecco run --genome <file> [options]
        gecco run --proteins <file> [options]

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
        input = self.args["--genome"] or self.args["--proteins"]
        if not os.path.exists(input):
            self.logger.error("could not locate input file: {!r}", input)
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
        self.logger.info("Running domain annotation")

        # Run PFAM HMM DB over ORFs to annotate with Pfam domains
        features = []
        for hmm in glob.glob(os.path.join(data.realpath("hmms"), "*.hmm.gz")):
            hmm_base = os.path.basename(hmm.replace(".hmm.gz", ""))
            self.logger.debug("Using HMM file {!r}", hmm_base)
            hmmer_out = os.path.join(out_dir, "hmmer", hmm_base)
            os.makedirs(hmmer_out, exist_ok=True)
            hmmer = HMMER(orf_file, hmmer_out, hmm, prodigal, self.args["--jobs"])
            features.append(hmmer.run().assign(hmm=hmm_base))

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


        # --- CRF ------------------------------------------------------------
        self.logger.info("Predicting cluster probabilities with the CRF model")
        crf = data.load("model/crf.model")

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
            found_clusters = refiner.find_clusters(
                subdf,
                method=self.args["--postproc"],
                prefix=sid,
            )
            if found_clusters:
                clusters.extend(found_clusters)

        if not clusters:
            self.logger.warning("No gene clusters were found")
            return 0

        # --- KNN ------------------------------------------------------------
        self.logger.info("Predicting BGC types")

        # Reformat training matrix
        self.logger.debug("Reading embedded training matrix")
        training_matrix = data.realpath("knn/domain_composition.tsv")
        train_df = pandas.read_csv(training_matrix, sep="\t", encoding="utf-8")
        train_comp = train_df.iloc[:,1:].values
        id_array = train_df["BGC_id"].values
        pfam_array = train_df.columns.values[1:]

        # Reformat type labels
        self.logger.debug("Reading embedded type labels")
        labels = data.realpath("knn/type_labels.tsv")
        types_df = pandas.read_csv(labels, sep="\t", encoding="utf-8")
        types_array = types_df["cluster_type"].values
        subtypes_array = types_df["subtype"].values

        # Calculate new domain composition
        self.logger.debug("Calulating domain composition for each cluster")
        new_comp = numpy.array([c.domain_composition(pfam_array) for c in clusters])

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
            csv.writer(f, dialect="excel-tab").writerow([
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
            ])
            for cluster, ty in zip(clusters, knn_pred):
                cluster.type, cluster.type_prob = ty
                cluster.write_to_file(f, long=True)

        # Write predicted cluster sequences to file
        for cluster in clusters:
            prots = []
            for p in SeqIO.parse(orf_file, "fasta"):
                if p.id in cluster.prot_ids:
                    p.description = f"{cluster.name} # {p.description}"
                    prots.append(p)

            prots_out = os.path.join(out_dir, f"{cluster.name}.proteins.faa")
            self.logger.debug("Writing proteins of {} to {!r}", cluster.name, prots_out)
            with open(prots_out, "w") as out:
                SeqIO.write(prots, out, "fasta")

        # Exit gracefully
        self.logger.info("Successfully found {} clusters!", len(clusters))
        return 0
