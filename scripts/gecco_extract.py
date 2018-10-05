import sys
import os
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
GECCO = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + "/..")
sys.path.append(GECCO)
import random
import argparse
import pickle
import pandas as pd
import numpy as np
from gecco.knn import ClusterKNN
from gecco.refine import ClusterRefiner
from gecco.interface import scripts_interface

SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])))
TRAINING_MATRIX = os.path.join(SCRIPT_DIR, "../data/knn/domain_composition.tsv")
LABELS = os.path.join(SCRIPT_DIR, "../data/knn/type_labels.tsv")

# MAIN
if __name__ == "__main__":
    args = scripts_interface()
    print(args)

    # Load data
    pfam_df = pd.read_csv(args.DATA, sep="\t", encoding="utf-8")
    try:
        pfam_df = pfam_df.sort_values(args.sort_cols).reset_index(drop=True)
    except KeyError as err:
        print("Colums could not be sorted.")

    # Extract clusters
    refiner = ClusterRefiner(threshold=args.thresh)
    clusters = []
    for sid, subdf in pfam_df.groupby(args.split_col):
        if len(subdf[args.group_col].unique()) < args.orfs:
            print(f"Skipping sequence {sid} because it is too short (< 5 ORFs).")
            continue
        found_clusters = refiner.find_clusters(
            subdf,
            method = args.post,
            prefix = sid
        )
        if found_clusters:
            clusters += found_clusters

    del pfam_df

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
    cluster_prots = np.hstack(np.array([c.prot_ids for c in clusters]))
    pfam_df["Y_pred"] = np.where(pfam_df[args.group_col].isin(cluster_prots), 1, 0)
    pfam_df.to_csv(f"{args.out}.refined.tsv", sep="\t", index=False, header=True)

    with open(f"{args.out}.clusters.tsv", "wt") as f:
        for c, t in zip(clusters, knn_pred):
            c.type = t[0]
            c.type_prob = t[1]
            c.write_to_file(f, long=True)
