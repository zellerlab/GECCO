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

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_predict.py /home/fleck/scripts/clust/test/test.embed.tsv /home/fleck/bin/orion/data/model/f5_eval_s_t50.crf.model  -o /home/fleck/scripts/clust/test/test


# FUNC
def interface():
    parser = argparse.ArgumentParser(description="Trains a model on a given set of samples and writes the model file and the formatted features to file.")

    parser.add_argument("DATA",
                        type=str,
                        metavar="<DATA>",
                        help="Pfam tables with training instances.")

    parser.add_argument("MODEL",
                        type=str,
                        metavar="<MODEL>",
                        help="Path to model.")

    parser.add_argument("-o", "--output-basename",
                        dest="out",
                        type=str,
                        default="CRF",
                        metavar="<basename>",
                        help="Basename for predictions.")

    parser.add_argument("-w", "--weight-col",
                        dest="w",
                        default="rev_i_Evalue",
                        type=str,
                        help="Column to be used as local weights on pfam domains.")

    parser.add_argument("--feature-type",
                        dest="feature_type",
                        type=str,
                        default="single",
                        help="How features should be extracted. 'Single', 'overlap' or on some grouping level ('group').")

    parser.add_argument("--overlap",
                        dest="overlap",
                        type=int,
                        default="2",
                        help="If overlapping features: How much overlap.")

    args = parser.parse_args()
    return args


# MAIN
if __name__ == "__main__":
    args = interface()

    data = args.DATA
    model_file = args.MODEL
    out_file = args.out

    C1 = 0.15
    C2 = 1.75

    weight_col = args.w
    feature_type = args.feature_type
    overlap = args.overlap

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    data_tbl = [s for _, s in data_tbl.groupby("BGC_id")]

    crf = ClusterCRF(data_tbl, "BGC",
        feature_cols = ["pfam"],
        weight_cols = [weight_col],
        feature_type = feature_type,
        overlap = overlap,
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2)

    with open(model_file, "rb") as f:
        model = pickle.load(f)
        crf.model = model

    pred_df = crf.predict_marginals()

    print(pred_df)

    pred_df.to_csv(out_file + ".pred.tsv", sep="\t", index=False, header=False)
