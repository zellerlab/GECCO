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
from scipy.stats import entropy
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

    data = "/Users/Jonas/Documents/Msc-Biotechnologie/masterarbeit-zeller/remote/scripts/clust/test/test.pred.tsv"
    thresh = 0.5
    split_col = "sequence_id"

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

    all_dom = list(cluster_list[0].domains) + list(cluster_list[1].domains) + list(cluster_list[2].domains)

    cluster_comp = np.array([c.domain_composition(all_possible=all_dom) for c in cluster_list])
