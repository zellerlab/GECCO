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

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_loto.py /home/fleck/scripts/clust/test/test.embed.tsv -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = crf_interface()

    data = args.DATA
    data_base = data.split("/")[-1].split(".")[0]
    out_file = args.out

    C1 = 0.15
    C2 = 1.75

    threads = args.threads
    if not threads:
        threads = multiprocessing.cpu_count()

    e_filter = args.e_filter
    truncate = args.truncate
    shuffle = args.shuffle
    weight_col = args.w
    feature_type = args.feature_type
    overlap = args.overlap

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    data_tbl = [s for _, s in data_tbl.groupby("BGC_id")]
    if shuffle:
        random.shuffle(data_tbl)

    crf = ClusterCRF(
        Y_col = "BGC",
        feature_cols = ["pfam"],
        weight_cols = [weight_col],
        feature_type = feature_type,
        overlap = overlap,
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2
    )

    results = crf.loto_cv(
        data_tbl,
        type_col = "BGC_type",
        threads = threads,
        e_filter = e_filter,
        truncate = truncate
    )

    result_df = (pd .concat(results)
                    .assign(c1 = C1,
                        c2 = C2,
                        feature_type = feature_type,
                        e_filter = e_filter,
                        overlap = overlap,
                        weight = weight_col,
                        truncate = truncate,
                        in_file = data_base,
                        cv_type = "LOTO")
                    .loc[ : , ["BGC", "BGC_id", "protein_id", "pfam", "pseudo_pos",
                        "p_pred", "c1", "c2", "feature_type", "e_filter", "overlap",
                        "weight", "truncate", "cv_type", "cv_round", "in_file"] ])

    # print(result_df)

    # Write results
    ext = "_loto" + ".pred.tsv"
    result_df.to_csv(out_file + ext, sep="\t", index=False, header=False)
