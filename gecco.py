#!/usr/bin/env python

#########################################################################################
#                                                                                       #
#                                           GECCO                                       #
#                                  GEne Cluster prediction                              #
#                               with COnditional random fields                          #
#                                                                                       #
#                                       MAIN SCRIPT                                     #
#                                                                                       #
#   Author: Jonas Simon Fleck (jonas.simon.fleck@gmail.com)                             #
#                                                                                       #
#########################################################################################

import os
import sys
import pickle
import logging
import argparse
import warnings
import subprocess
import multiprocessing
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import numpy as np
import pandas as pd
from gecco.hmmer import HMMER
from gecco.orf import ORFFinder
from gecco.crf import ClusterCRF
from gecco.knn import ClusterKNN
from gecco.refine import ClusterRefiner
from gecco.interface import main_interface

# CONST
SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])))
PFAM = open(os.path.join(SCRIPT_DIR, "data/db_config.txt")).readlines()[0].strip()
MODEL = os.path.join(SCRIPT_DIR, "data/model/feat_v8_param_v2.crf.model")
TRAINING_MATRIX = os.path.join(SCRIPT_DIR, "data/knn/domain_composition.tsv")
LABELS = os.path.join(SCRIPT_DIR, "data/knn/type_labels.tsv")

# MAIN
if __name__ == "__main__":
    # PARAMS
    args = main_interface()

    # Make out directory
    out_dir = args.out
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Set up logging
    logging.basicConfig(
        level = logging.INFO,
        format = "%(asctime)s [%(levelname)s]  %(message)s",
        handlers = [
            logging.FileHandler(os.path.join(out_dir, "gecco.log")),
            logging.StreamHandler()
    ])

    logging.info(f"GECCO is running with these parameters:\n{args.__dict__}")

    genome = args.GENOME
    base = ".".join(os.path.basename(genome).split(".")[:-1])

    e_filter = min(1, args.e_filter)
    threads = args.threads
    if not threads:
        threads = multiprocessing.cpu_count()


    # PRODIGAL
    logging.info("Predicting ORFs with PRODIGAL.")

    prodigal_out = os.path.join(out_dir, "prodigal/")
    if not os.path.exists(prodigal_out):
        os.makedirs(prodigal_out)

    # Extract ORFs from genome
    prodigal = ORFFinder(genome, prodigal_out, method="prodigal")
    orf_file = prodigal.run()


    # HMMER
    logging.info("Running Pfam domain annotation.")

    hmmer_out = os.path.join(out_dir, "hmmer/")
    if not os.path.exists(hmmer_out):
        os.makedirs(hmmer_out)

    # Run PFAM HMM DB over ORFs to annotate with Pfam domains
    hmmer = HMMER(orf_file, hmmer_out, hmms=PFAM)
    pfam_df = hmmer.run()

    # Filter i-Evalue
    pfam_df = pfam_df[pfam_df["i_Evalue"] < e_filter]
    # Reformat pfam IDs
    pfam_df = pfam_df.assign(
        pfam = pfam_df["pfam"].str.replace(r"(PF\d+)\.\d+", lambda m: m.group(1))
    )

    # Write feature table to file
    feat_out = os.path.join(out_dir, base + ".features.tsv")
    pfam_df.to_csv(feat_out, sep="\t", index=False)


    # CRF
    logging.info("Prediction of cluster probabilities with the CRF model.")

    # Load model from file
    with open(MODEL, "rb") as f:
        crf = pickle.load(f)

    # Split input dataframe and predict marginal probabilitites
    pfam_df = [seq for _, seq in pfam_df.groupby("sequence_id")]
    pfam_df = crf.predict_marginals(data=pfam_df)

    # Write predictions to file
    pred_out = os.path.join(out_dir, base + ".pred.tsv")
    pfam_df.to_csv(pred_out, sep="\t", index=False)


    # REFINE
    logging.info("Extracting clusters.")

    refiner = ClusterRefiner(threshold=args.thresh)

    clusters = []
    for sid, subdf in pfam_df.groupby("sequence_id"):
        if len(subdf["protein_id"].unique()) < 5:
            logging.warning(
                f"Skipping sequence {sid} because it is too short (< 5 ORFs).")
            continue
        found_clusters = refiner.find_clusters(
            subdf,
            method = args.post,
            prefix = sid
        )
        if found_clusters:
            clusters += found_clusters

    del pfam_df

    if not clusters:
        logging.warning("Unfortunately, no clusters were found. Exiting now.")
        sys.exit()


    # KNN
    logging.info("Prediction of BGC types.")

    # Reformat training matrix
    train_df = pd.read_csv(TRAINING_MATRIX, sep="\t", encoding="utf-8")
    train_comp = train_df.iloc[:,1:].values
    id_array = train_df["BGC_id"].values
    pfam_array = train_df.columns.values[1:]

    # Reformant type labels
    types_df = pd.read_csv(LABELS, sep="\t", encoding="utf-8")
    types_array = types_df["cluster_type"].values
    subtypes_array = types_df["subtype"].values

    # Calculate domain composition for all new found clusters
    new_comp = np.array(
        [c.domain_composition(all_possible=pfam_array) for c in clusters]
    )

    # Inititate kNN and predict types
    knn = ClusterKNN(metric=args.dist, n_neighbors=args.k)
    knn_pred = knn.fit_predict(train_comp, new_comp, y=types_array)

    # Write predicted clusters to file
    cluster_out = os.path.join(out_dir, base + ".clusters.tsv")
    with open(cluster_out, "wt") as f:
        for c, t in zip(clusters, knn_pred):
            c.type = t[0]
            c.type_prob = t[1]
            c.write_to_file(f, long=True)

    logging.info("DONE.\n")
