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
from orion.interface import refine_interface

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_refine.py /home/fleck/scripts/clust/test/test.pred.tsv -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = refine_interface()

    data = args.DATA
    out_file = args.out
    thresh = args.thresh
    split_col = args.split_col

    print(args)

    data_df = pd.read_csv(data, sep="\t", encoding="utf-8")

    refiner = ClusterRefiner(
        threshold = thresh
    )

    cluster_list = []
    for sid, df in data_df.groupby(split_col):
        clusters = refiner.find_clusters(
            df,
            method = "antismash",
            prefix = sid
        )
        if clusters:
            cluster_list += clusters

    with open(out_file + ".clusters.tsv", "wt") as f:
        for c in clusters:
            c.write_to_file(f)
