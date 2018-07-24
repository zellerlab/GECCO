import sys
import os
ORION = os.path.abspath(os.path.dirname(sys.argv[0]) + "/..")
sys.path.append(ORION)
import random
import argparse
import pickle
import pandas as pd
import numpy as np
import multiprocessing
from itertools import product
from orion.crf import ClusterCRF
from orion.interface import crf_interface

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_train.py /home/fleck/scripts/clust/test/test.embed.tsv -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = crf_interface()

    data = args.DATA
    out_file = args.out

    C1 = 0.15
    C2 = 1.75

    trunc = args.truncate
    shuffle = args.shuffle
    weight_col = args.w
    feature_type = args.feature_type
    overlap = args.overlap

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    data_tbl = [s for _, s in data_tbl.groupby("BGC_id")]
    if shuffle:
        random.shuffle(data_tbl)

    crf = ClusterCRF(data_tbl, "BGC",
        feature_cols = ["pfam"],
        weight_cols = [weight_col],
        feature_type = feature_type,
        overlap = overlap,
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2)

    if trunc:
        crf.truncate(trunc)

    crf.fit()

    with open(out_file + ".crf.model", "wb") as f:
        pickle.dump(crf.model, f, protocol=2)

    with open(out_file + ".trans.tsv", "wt") as f:
        f.write("from\tto\tweight\n")
        for (label_from, label_to), weight in crf.model.transition_features_.items():
            f.write("{0}\t{1}\t{2}\n".format(label_from, label_to, weight))

    with open(out_file + ".state.tsv", "wt") as f:
        f.write("attr\tlabel\tweight\n")
        for (attr, label), weight in crf.model.state_features_.items():
            f.write("{0}\t{1}\t{2}\n".format(attr, label, weight))
