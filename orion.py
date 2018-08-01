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
import argparse
import pickle
import subprocess
import pandas as pd
import numpy as np
from orion.utils import convert_hmmer
from orion.crf import ClusterCRF
from orion.interface import main_interface
from orion.preprocessing import compute_features

# CONST
SCRIPT_DIR = os.path.abspath(os.path.dirname(sys.argv[0]))
PFAM = open(os.path.join(SCRIPT_DIR, "data/db_config.txt")).readlines()[0].strip()
MODEL = os.path.join(SCRIPT_DIR, "data/model/f5_eval_p_t50.crf.model")
C1 = 0.15
C2 = 1.75

# MAIN
if __name__ == "__main__":

    # PARAMS
    args = main_interface()

    print("Running ORION with these parameters:")
    print(args)

    fasta = args.FASTA
    base = os.path.basename(fasta).split(".")[0]

    out_dir = args.out
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    e_filter = max(1, args.e_filter)
    weight_type = args.weight_type

    # PRODIGAL

    print("Running ORF prediction using PRODIGAL...")

    prodigal_out = os.path.join(out_dir, "prodigal/")
    if not os.path.exists(prodigal_out):
        os.makedirs(prodigal_out)

    coords_out = os.path.join(prodigal_out, base + ".coords.sco")
    genes_out = os.path.join(prodigal_out, base + ".genes.fna")
    proteins_out = os.path.join(prodigal_out, base + ".proteins.faa")
    log_out = os.path.join(prodigal_out, base + ".prodigal.log")
    subprocess.run(["prodigal",
        "-i", fasta,
        "-o", coords_out,
        "-f", "sco",
        "-d", genes_out,
        "-a", proteins_out,
        "-p", "meta"],
        stdout = open(log_out, "wt"),
        stderr = open(log_out, "wt"))

    # HMMER

    print("Running Pfam domain annotation...")

    with open(proteins_out, "rt") as f:
        gene_order = [l.split()[0][1:] for l in f if l.startswith(">")]

    hmmer_out = os.path.join(out_dir, "hmmer/")
    if not os.path.exists(hmmer_out):
        os.makedirs(hmmer_out)

    dom_out = os.path.join(hmmer_out, base + ".hmmer.dom")
    log_out = os.path.join(hmmer_out, base + ".hmmer.log")
    subprocess.run(["hmmsearch",
        "--domtblout", dom_out, PFAM, proteins_out],
        stdout = open(log_out, "wt"),
        stderr = open(log_out, "wt"))

    tsv_out = os.path.join(hmmer_out, base + ".hmmer.tsv")
    convert_hmmer(dom_out, tsv_out)

    # Format feature table
    pfam_df = pd.read_csv(tsv_out, sep = "\t")
    pfam_df["protein_id"] = pd.Categorical(pfam_df["protein_id"], gene_order)
    pfam_df = pfam_df[pfam_df["i_Evalue"] < e_filter]
    pfam_df = compute_features(pfam_df, weight_type=weight_type)

    # Write feature table to file
    feat_out = os.path.join(out_dir, base + ".features.tsv")
    pfam_df.to_csv(feat_out, sep="\t", index=False)


    # Predict

    print("Running Cluster prediction...")

    crf = ClusterCRF(
        data = [pfam_df],
        group_col = "protein_id",
        feature_type = "group",
        feature_cols = ["pfam"],
        weight_cols = ["noweight"],
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2
    )

    with open(MODEL, "rb") as f:
        model = pickle.load(f)
        crf.model = model

    pfam_df = crf.predict_marginals()

    # Write predictions to file
    pred_out = os.path.join(out_dir, base + ".pred.tsv")
    pfam_df.to_csv(pred_out, sep="\t", index=False)

    print("DONE.")
