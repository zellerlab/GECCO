import sys
import os
MAIN = os.path.abspath(os.path.dirname(sys.argv[0]) + "/..") + "/"
sys.path.append(MAIN)
import random
import argparse
import pickle
import pandas as pd
import numpy as np
import multiprocessing
from itertools import product
from main.crf import ClusterCRF

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_train.py /home/fleck/scripts/clust/test/test.embed.tsv -o /home/fleck/scripts/clust/test/test


# FUNC
def interface():
    parser = argparse.ArgumentParser(description="Trains a model on a given set of samples and writes the model file and the formatted features to file.")

    parser.add_argument("DATA",
                        type=str,
                        metavar="<DATA>",
                        help="Pfam tables with training instances.")

    parser.add_argument("-o", "--output-basename",
                        dest="out",
                        type=str,
                        default="CRF",
                        metavar="<basename>",
                        help="Basename for model binary and tsv.")

    parser.add_argument("-w", "--weight-col",
                        dest="w",
                        default="pseudo_norm",
                        type=str,
                        help="Column to be used as local weights on pfam domains.")

    parser.add_argument("--feature-type",
                        dest="feature_type",
                        type=str,
                        default="group",
                        help="How features should be extracted. 'Single', 'overlap' or on some grouping level ('group').")

    parser.add_argument("--truncate",
                        dest="truncate",
                        type=int,
                        help="Training set will be truncated to this length.")

    parser.add_argument("--overlap",
                        dest="overlap",
                        type=int,
                        default="2",
                        help="If overlapping features: How much overlap.")

    parser.add_argument("--no-shuffle",
                        dest="shuffle",
                        action="store_false",
                        help="Switch to turn of shuffling of the data before doing CV.")

    args = parser.parse_args()
    return args


# MAIN
if __name__ == "__main__":
    args = interface()

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
