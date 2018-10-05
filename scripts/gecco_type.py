import os
import sys
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
GECCO = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + "/..")
sys.path.append(GECCO)
import random
import argparse
import pickle
import os
import pandas as pd
import numpy as np
from gecco.knn import ClusterKNN
from gecco.bgc import BGC
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
    clusters = list()
    with open(args.DATA, "rt") as f:
        for line in f:
            line = line.split()
            cluster = BGC(
                seq_id = line[1],
                domains = line[3].split(","),
                weights = map(float, line[4].split(",")),
                name = os.path.abspath(args.DATA)
            )
            clusters.append(cluster)

    # Reformat training matrix
    train_df = pd.read_csv(TRAINING_MATRIX, sep="\t", encoding="utf-8")
    train_comp = train_df.iloc[:,1:].values
    id_array = train_df["BGC_id"].values
    pfam_array = train_df.columns.values[1:]

    # Reformant type labels
    types_df = pd.read_csv(LABELS, sep="\t", encoding="utf-8")
    types_array = types_df["cluster_type"].values
    subtypes_array = types_df["subtype"].values

    # Calculate domain composition for all new found cluster
    new_comp = np.array(
        [c.domain_composition(all_possible=pfam_array) for c in clusters]
    )

    # Inititate kNN and predict types
    knn = ClusterKNN(metric=args.dist, n_neighbors=args.k)
    knn_pred = knn.fit_predict(train_comp, new_comp, y=types_array)

    # Write predicted clusters to file
    with open(f"{args.out}.clusters.tsv", "wt") as f:
        for c, t in zip(clusters, knn_pred):
            c.type = t[0]
            c.type_prob = t[1]
            c.write_to_file(f)
