#!/usr/bin/env python

#########################################################################################
#                                                                                       #
#                                           ORION                                       #
#                           predicting biosynthetic gene clusters                       #
#                              using conditional random fields                          #
#                                                                                       #
#                                       MAIN SCRIPT                                     #
#                                                                                       #
#   Author: Jonas Simon Fleck (jonas.simon.fleck@gmail.com)                             #
#                                                                                       #
#########################################################################################

import os
import sys
import pickle
import argparse
import subprocess
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import pandas as pd
import numpy as np
from orion.orf import ORFFinder
from orion.hmmer import HMMER
from orion.crf import ClusterCRF
from orion.knn import ClusterKNN
from orion.refine import ClusterRefiner
from orion.interface import main_interface
from orion.preprocessing import compute_features

# CONST
SCRIPT_DIR = os.path.abspath(os.path.dirname(sys.argv[0]))
PFAM = open(os.path.join(SCRIPT_DIR, "data/db_config.txt")).readlines()[0].strip()
MODEL = os.path.join(SCRIPT_DIR, "data/model/f5_eval_p_t50.crf.model")
TRAINING_MATRIX = os.path.join(SCRIPT_DIR, "data/knn/domain_composition.tsv")
LABELS = os.path.join(SCRIPT_DIR, "data/knn/type_labels.tsv")

# MAIN
if __name__ == "__main__":

    # PARAMS
    args = main_interface()

    sys.stdout.write("Running ORION with these parameters:" + "\n")
    sys.stdout.write(str(args) + "\n")

    fasta = args.FASTA
    base = os.path.basename(fasta).split(".")[0]

    out_dir = args.out
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    e_filter = max(1, args.e_filter)
    weight_type = args.weight_type


    # PRODIGAL

    sys.stdout.write("Running ORF prediction using PRODIGAL..." + "\n")

    prodigal_out = os.path.join(out_dir, "prodigal/")
    if not os.path.exists(prodigal_out):
        os.makedirs(prodigal_out)

    prodigal = ORFFinder(fasta, prodigal_out, method="prodigal")
    orf_file = prodigal.run()


    # HMMER

    sys.stdout.write("Running Pfam domain annotation..." + "\n")

    hmmer_out = os.path.join(out_dir, "hmmer/")
    if not os.path.exists(hmmer_out):
        os.makedirs(hmmer_out)

    hmmer = HMMER(orf_file, hmmer_out, hmms=PFAM)
    pfam_df = hmmer.run()

    # Format feature table
    pfam_df = pfam_df[pfam_df["i_Evalue"] < e_filter]
    # pfam_df = compute_features(pfam_df, weight_type=weight_type)

    # Write feature table to file
    feat_out = os.path.join(out_dir, base + ".features.tsv")
    pfam_df.to_csv(feat_out, sep="\t", index=False)


    # CRF

    sys.stdout.write("Running cluster prediction..." + "\n")

    with open(MODEL, "rb") as f:
        crf = pickle.load(f)

    ### TEMPORARY HACK, HAVE TO REPLACE MODEL ###
    crf.weights = [1]
    #############################################

    pfam_df = crf.predict_marginals(data=[pfam_df])

    # Write predictions to file
    pred_out = os.path.join(out_dir, base + ".pred.tsv")
    pfam_df.to_csv(pred_out, sep="\t", index=False)


    # REFINE

    sys.stdout.write("Extracting and refining clusters..." + "\n")

    refiner = ClusterRefiner(threshold=0.5)
    clusters = refiner.find_clusters(
        pfam_df,
        method = "antismash",
        prefix = sid
    )

    if not clusters:
        sys.stdout.write("Unfortunately, no clusters were found. Exiting now.")
        sys.exit()


    # KNN

    sys.stdout.write("Running cluster type prediction..." + "\n")

    train_df = pd.read_csv(TRAINING_MATRIX, sep="\t", encoding="utf-8")
    train_comp = train_df.iloc[:,1:].values
    id_array = train_df["BGC_id"].values
    pfam_array = train_df.columns.values[1:]

    types_df = pd.read_csv(LABELS, sep="\t", encoding="utf-8")
    types_array = types_df["cluster_type"].values
    subtypes_array = types_df["subtype"].values

    new_comp = np.array(
        [c.domain_composition(all_possible=pfam_array) for c in clusters]
    )

    knn = ClusterKNN(metric="jsd", n_neighbors=1)
    knn_pred = knn.fit_predict(train_comp, new_comp, y=types_array)

    cluster_out = os.path.join(out_dir, base + ".clusters.tsv")
    with open(cluster_out, "wt") as f:
        for c, t in zip(clusters, knn_pred):
            c.type = t
            c.write_to_file(f)

    sys.stdout.write("DONE." + "\n")
