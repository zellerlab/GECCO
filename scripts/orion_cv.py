import sys
import os
ORION = os.path.abspath(os.path.dirname(sys.argv[0]) + "/..")
sys.path.append(ORION)
import random
import argparse
import pandas as pd
import numpy as np
import multiprocessing
from itertools import product
from orion.crf import ClusterCRF
from orion.interface import crf_interface
from orion.utils import coerce_numeric

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_cv.py /home/fleck/scripts/clust/test/test.embed.tsv -o /home/fleck/scripts/clust/test/test -t1 --split-col BGC_id --sort-col BGC_id start --folds 2

# MAIN
if __name__ == "__main__":
    args = crf_interface()

    data = args.DATA
    data_base = data.split("/")[-1].split(".")[0]
    out_file = args.out


    threads = args.threads
    if not threads:
        threads = multiprocessing.cpu_count()

    C1 = args.C1
    C2 = args.C1
    weight_col = [coerce_numeric(w) for w in args.w]
    y_col = args.y
    feature_col = args.feat
    group_col = args.group_col
    strat_col = args.strat_col
    split_col = args.split_col
    trunc = args.truncate
    splits = args.splits
    shuffle = args.shuffle
    feature_type = args.feature_type
    overlap = args.overlap
    e_filter = args.e_filter

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    data_tbl = [s for _, s in data_tbl.groupby(split_col)]
    if shuffle:
        random.shuffle(data_tbl)

    crf = ClusterCRF(
        Y_col = y_col,
        feature_cols = feature_col,
        weight_cols = weight_col,
        feature_type = feature_type,
        overlap = overlap,
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2
    )

    results = crf.cv(
        data_tbl,
        k = splits,
        threads = threads,
        e_filter = e_filter,
        trunc = trunc,
        strat_col = strat_col
    )

    result_df = (pd .concat(results)
                    .assign(c1 = C1,
                        c2 = C2,
                        feature_type = feature_type,
                        e_filter = e_filter,
                        overlap = overlap,
                        weight = ",".join(map(str, weight_col)),
                        feature = ",".join(feature_col),
                        truncate = trunc,
                        in_file = data_base,
                        cv_round = "all",
                        cv_type = "10-fold")
                    .loc[ : , ["BGC", "BGC_id", "protein_id", "pfam", "pseudo_pos",
                        "p_pred", "c1", "c2", "feature_type", "e_filter", "overlap",
                        "weight", "truncate", "cv_type", "cv_round", "in_file"] ])

    # print(result_df)

    # Write results
    ext = "_cv" + ".pred.tsv"
    result_df.to_csv(out_file + ext, sep="\t", index=False, header=False)
