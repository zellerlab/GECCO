# coding: utf-8

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

    summary = "predict Biosynthetic Gene Clusters from a genome file."
    doc = f"""
    gecco run - {summary}

    Usage:
        gecco run --genome <file>  [options]
        gecco run --proteins <file> [options]
        gecco run (-h | --help)

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
        -e <e>, --e-filter <e>        the e-value cutoff for PFam domains to
                                      be included [default: 1e-5]
        -m <m>, --threshold <m>       the probability threshold for cluster
                                      detection. Default depends on the
                                      post-processing method (0.4 for gecco,
                                      0.6 for antismash).
        -k <n>, --neighbors <n>       the number of neighbors to use for
                                      kNN type prediction [default: 5]
        -d <d>, --distance <d>        the distance metric to use for kNN type
                                      prediction. [default: jensenshannon]
        --postproc <method>           the method to use for cluster extraction
                                      (antismash or gecco). [default: gecco]
    """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check value of numeric arguments
        self.args["--neighbors"] = int(self.args["--neighbors"])
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
            self.logger.critical("could not locate input file: {!r}", input)
            return 1

        return None

    def __call__(self) -> int:
        # Check CLI arguments
        retcode = self._check()
        if retcode is not None:
            return retcode

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
            orf_file = self.args["--protein"]
            base, _ = os.path.splitext(os.path.basename(orf_file))
            prodigal = False

        # --- HMMER ----------------------------------------------------------
        self.logger.info("Running PFam domain annotation")
        hmmer_out = os.path.join(out_dir, "hmmer")
        os.makedirs(hmmer_out, exist_ok=True)

        # Run PFAM HMM DB over ORFs to annotate with Pfam domains
        hmms = data.realpath("hmms/Pfam-A.hmm.gz")
        hmmer = HMMER(orf_file, hmmer_out, hmms, prodigal, self.args["--jobs"])
        pfam_df = hmmer.run()

        # Filter i-evalue
        self.logger.debug("Filtering results with e-value under {}", self.args["--e-filter"])
        pfam_df = pfam_df[pfam_df["i_Evalue"] < self.args["--e-filter"]]

        # Reformat pfam IDs
        pfam_df = pfam_df.assign(
            pfam=pfam_df["pfam"].str.replace(r"(PF\d+)\.\d+", lambda m: m.group(1))
        )

        # Write feature table to file
        feat_out = os.path.join(out_dir, f"{base}.features.tsv")
        self.logger.debug("Writing feature table to {!r}", feat_out)
        pfam_df.to_csv(feat_out, sep="\t", index=False)


        # --- CRF ------------------------------------------------------------
        self.logger.info("Predicting cluster probabilities with the CRF model")
        with data.open("model/feat_v8_param_v2.crf.model", "rb") as bin:
            crf = pickle.load(bin)

        # If extracted from genome, split input dataframe into sequence
        if prodigal:
            pfam_df = [seq for _, seq in pfam_df.groupby("sequence_id")]
        else:
            pfam_df = [pfam_df]
        pfam_df = crf.predict_marginals(data=pfam_df)

        # Write predictions to file
        pred_out = os.path.join(out_dir, f"{base}.pred.tsv")
        pfam_df.to_csv(pred_out, sep="\t", index=False)


        # --- REFINE ---------------------------------------------------------
        self.logger.info("Extracting clusters")
        refiner = ClusterRefiner(threshold=self.args["--threshold"])

        clusters = []
        for sid, subdf in pfam_df.groupby("sequence_id"):
            if len(subdf["protein_id"].unique()) < 5:
                self.logger.warn("Skipping sequence {} because it is too short", sid)
                continue
            found_clusters = refiner.find_clusters(
                subdf,
                method=self.args["--postproc"],
                prefix=sid,
            )
            if found_clusters:
                clusters.extend(found_clusters)

        if not clusters:
            self.logger.warning("Unfortunately, no clusters were found.")
            return 0

        # --- KNN ------------------------------------------------------------
        self.logger.info("Predicting BGC types")

        # Reformat training matrix
        training_matrix = data.realpath("knn/domain_composition.tsv")
        train_df = pandas.read_csv(training_matrix, sep="\t", encoding="utf-8")
        train_comp = train_df.iloc[:,1:].values
        id_array = train_df["BGC_id"].values
        pfam_array = train_df.columns.values[1:]

        # Reformat type labels
        labels = data.realpath("knn/type_labels.tsv")
        types_df = pandas.read_csv(labels, sep="\t", encoding="utf-8")
        types_array = types_df["cluster_type"].values
        subtypes_array = types_df["subtype"].values

        # Calculate new domain composition
        new_comp = numpy.array(
            [c.domain_composition(all_possible=pfam_array) for c in clusters]
        )

        # Inititate kNN and predict types
        knn = ClusterKNN(metric=self.args["--distance"], n_neighbors=self.args["--neighbors"])
        knn_pred = knn.fit_predict(train_comp, new_comp, y=types_array)

        # --- RESULTS --------------------------------------------------------
        self.logger.info("Writing final results file")

        # Write predicted cluster coordinates to file
        cluster_out = os.path.join(out_dir, f"{base}.clusters.tsv")
        with open(cluster_out, "wt") as f:
            for cluster, ty in zip(clusters, knn_pred):
                cluster.type = ty[0]
                cluster.type_prob = ty[1]
                cluster.write_to_file(f, long=True)

        # Write predicted cluster sequences to file
        for cluster in clusters:
            prots = cluster.prot_ids
            cid = cluster.name
            prot_list = []
            proteins = SeqIO.parse(orf_file, "fasta")
            for p in proteins:
                if p.id in prots:
                    p.description = f"{cid} # {p.description}"
                    prot_list.append(p)
            with open(os.path.join(out_dir, f"{cid}.proteins.faa"), "wt") as out:
                SeqIO.write(prot_list, out, "fasta")

        # Exit gracefully
        self.logger.info("Successfully found {} clusters!", len(clusters))
        return 0
