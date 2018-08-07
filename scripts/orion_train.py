import sys
import os
ORION = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + "/..")
sys.path.append(ORION)
import random
import argparse
import pickle
import pandas as pd
import numpy as np
import multiprocessing
from itertools import product
from orion.crf import ClusterCRF
from orion.interface import scripts_interface
from orion.preprocessing import truncate
from orion.utils import coerce_numeric

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_train.py /home/fleck/scripts/clust/test/test.embed.tsv -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = scripts_interface()

    data = args.DATA
    out_file = args.out

    C1 = args.C1
    C2 = args.C1
    weight_col = [coerce_numeric(w) for w in args.w]
    y_col = args.y
    feature_col = args.feat
    group_col = args.group_col
    split_col = args.split_col
    trunc = args.truncate
    shuffle = args.shuffle
    feature_type = args.feature_type
    overlap = args.overlap

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    data_tbl = [s for _, s in data_tbl.groupby(split_col)]
    if shuffle:
        random.shuffle(data_tbl)

    crf = ClusterCRF(
        Y_col = y_col,
        feature_cols = feature_col,
        group_col = group_col,
        weight_cols = weight_col,
        feature_type = feature_type,
        overlap = overlap,
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2
    )

    if trunc:
        data_tbl = [truncate(df, trunc, Y_col=y_col, grouping=group_col)
            for df in data_tbl]

    crf.fit(data=data_tbl)

    with open(out_file + ".crf.model", "wb") as f:
        pickle.dump(crf, f, protocol=2)

    with open(out_file + ".trans.tsv", "wt") as f:
        f.write("from\tto\tweight\n")
        for (label_from, label_to), weight in crf.model.transition_features_.items():
            f.write("{0}\t{1}\t{2}\n".format(label_from, label_to, weight))

    with open(out_file + ".state.tsv", "wt") as f:
        f.write("attr\tlabel\tweight\n")
        for (attr, label), weight in crf.model.state_features_.items():
            f.write("{0}\t{1}\t{2}\n".format(attr, label, weight))
