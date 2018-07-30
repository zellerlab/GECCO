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
SCRIPT_DIR = os.path.abspath(os.path.dirname(sys.argv[0])) + "/"
PFAM = os.environ["PFAM"]
MODEL = SCRIPT_DIR + "data/model/f5_eval_p_t50.crf.model"
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
    if not out_dir.endswith("/"):
        out_dir += "/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    e_filter = max(1, args.e_filter)
    weight_type = args.weight_type

    # PRODIGAL

    print("Running ORF prediction using PRODIGAL...")

    prodigal_out = out_dir + "prodigal/"
    if not os.path.exists(prodigal_out):
        os.mkdir(prodigal_out)

    coords_out = prodigal_out + base + ".coords.sco"
    genes_out = prodigal_out + base + ".genes.fna"
    proteins_out = prodigal_out + base + ".proteins.faa"
    log_out = prodigal_out + ".prodigal.log"
    subprocess.run(["prodigal",
        "-i", fasta,
        "-o", coords_out,
        "-f", "sco",
        "-d", genes_out,
        "-a", proteins_out,
        "-p", "meta"],
        stdout = open(log_out, "wt"),
        stderr = open(log_out, "wt"))

    print("DONE")

    # HMMER

    print("Running Pfam domain annotation...")

    with open(proteins_out, "rt") as f:
        gene_order = [l.split()[0][1:] for l in f if l.startswith(">")]

    print(gene_order)

    hmmer_out = out_dir + "hmmer/"
    if not os.path.exists(hmmer_out):
        os.mkdir(hmmer_out)

    dom_out = hmmer_out + base + ".hmmer.dom"
    log_out = hmmer_out + ".hmmer.log"
    subprocess.run(["hmmsearch",
        "--domtblout", dom_out, PFAM, proteins_out],
        stdout = open(log_out, "wt"),
        stderr = open(log_out, "wt"))

    tsv_out = hmmer_out + base + ".hmmer.tsv"
    convert_hmmer(dom_out, tsv_out)

    print("DONE")

    # Format feature table
    pfam_df = pd.read_csv(tsv_out, sep = "\t")
    pfam_df["protein_id"] = pd.Categorical(pfam_df["protein_id"], gene_order)
    pfam_df = pfam_df[pfam_df["i_Evalue"] < e_filter]
    pfam_df = compute_features(pfam_df, weight_type=weight_type)

    # Write feature table to file
    feat_out = out_dir + base + ".features.tsv"
    pfam_df.to_csv(feat_out, sep="\t", index=False)


    # Predict

    print("Running Cluster prediction.")

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

    pred_df = crf.predict_marginals()

    # Write predictions to file
    pred_df = out_dir + base + ".pred.tsv"
    pfam_df.to_csv(pred_df, sep="\t", index=False)
