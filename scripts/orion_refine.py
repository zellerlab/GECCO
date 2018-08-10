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
from orion.refine import ClusterRefiner
from orion.interface import scripts_interface

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_refine.py /home/fleck/scripts/clust/test/test.pred.tsv -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = scripts_interface()

    data = args.DATA
    out_file = args.out

    print(args)

    data_df = pd.read_csv(data, sep="\t", encoding="utf-8")

    refiner = ClusterRefiner(
        threshold = args.thresh
    )

    cluster_list = []
    if args.split_col:
        for sid, df in data_df.groupby(args.split_col):
            clusters = refiner.find_clusters(
                df,
                method = args.post,
                prefix = sid
            )
            if clusters:
                cluster_list += clusters
    else:
        cluster_list = refiner.find_clusters(
            df,
            method = args.post,
            prefix = sid
        )

    cluster_prots = np.hstack(np.array([c.prot_ids for c in cluster_list]))
    data_df["Y_pred"] = np.where(data_df[args.group_col].isin(cluster_prots), 1, 0)

    data_df.to_csv(out_file + ".refined.tsv", sep="\t", index=False, header=True)

    with open(out_file + ".clusters.tsv", "wt") as f:
        for c in cluster_list:
            c.write_to_file(f, short=True)
