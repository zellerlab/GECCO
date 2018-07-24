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
import subprocess
import pandas as pd
import numpy as np
from main.utils import convert_hmmer
from main.crf import ClusterCRF

# CONST
SCRIPT_DIR = os.path.abspath(os.path.dirname(sys.argv[0])) + "/"
PFAM = SCRIPT_DIR + "data/pfam/Pfam-A.hmm"
MODEL = SCRIPT_DIR + "data/model/"

# FUNC
def interface():
    parser = argparse.ArgumentParser(description="Predicts biosynthetic gene clusters from a sorted FASTA file.")

    parser.add_argument("FASTA",
                        type=str,
                        metavar="<FASTA>",
                        help="FASTA file with proteins.")

    parser.add_argument("-o", "--output-dir",
                        dest="out",
                        type=str,
                        default="./",
                        metavar="<out_dir>",
                        help="Output directory.")

    parser.add_argument("--e-filter",
                        dest="e_filter",
                        type=float,
                        default="1e-5",
                        metavar="<e_filter>",
                        help="E-value cutoff for pfam domains to be included.")

    args = parser.parse_args()
    return args


# MAIN
if __name__ == "__main__":

    print(SCRIPT_DIR)

    args = interface()

    fasta = args.FASTA
    base = os.path.basename(fasta).split(".")[0]

    out_dir = args.out
    if not out_dir.endswith("/"):
        out_dir += "/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    e_filter = max(1, args.e_filter)

    with open(fasta, "rt") as f:
        gene_order = [l.split()[0][1:] for l in f if l.startswith(">")]

    print(gene_order)

    # HMMER


    hmmer_out = out_dir + "hmmer/"
    if not os.path.exists(hmmer_out):
        os.mkdir(hmmer_out)

    dom_out = hmmer_out + base + ".hmmer.dom"
    log_out = hmmer_out + "hmmer.log"
    subprocess.run(["hmmsearch", "--domtblout", dom_out, PFAM, fasta],
        stdout = open(log_out, "wt"))

    tsv_out = hmmer_out + base + ".hmmer.tsv"
    convert_hmmer(dom_out, tsv_out)

    # Format feature table
    pfam_df = pd.read_csv(tsv_out, sep = "\t")
    pfam_df["protein_id"] = pd.Categorical(pfam_df["protein_id"], gene_order)
    pfam_df = pfam_df[pfam_df["i_Evalue"] < e_filter]
    pfam_df = (pfam_df
        .sort_values(by = "protein_id")
        .assign(
            pseudo_pos = range(len(pfam_df)),
            rev_Evalue = 1 - pfam_df["i_Evalue"],
            log_Evalue = [min(- np.log10(n), 300) for n in pfam_df["i_Evalue"]]
            )
        .groupby("pfam", sort = False)
        .apply(lambda x: x.assign(
            supnorm_Evalue = (x["log_Evalue"] / x["log_Evalue"].max())
        ))
        .reset_index(drop = True)
        .groupby("protein_id", sort = False)
        .apply(lambda x: x.assign(
            supnorm_Evalue_prot = x["supnorm_Evalue"] / x["supnorm_Evalue"].sum()
        ))
        .reset_index(drop = True)
        .sort_values("pseudo_pos")
    )

    # Write feature table to file
    feat_out = out_dir + base + ".features.tsv"
    pfam_df.to_csv(feat_out, sep="\t", index=False)

    print(pfam_df)

    # Predict
    crf = ClusterCRF(
        data = pfam_df,
        model_filename = MODEL,
        group_col = "protein_id",
        feature_type = "group",
        feature_cols = ["pfam"],
        weight_cols = ["supnorm_Evalue_prot"]
    )
    crf.predict()
    marginal_probs = self.crf.predict_marginals()
