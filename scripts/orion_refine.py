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
from orion.refine import ClusterRefiner
from orion.interface import refine_interface

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_predict.py /home/fleck/scripts/clust/test/test.embed.tsv /home/fleck/bin/orion/data/model/f5_eval_s_t50.crf.model  -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = refine_interface()

    data = args.DATA
    out_file = args.out
    thresh = args.thresh

    print(args)

    data_df = pd.read_csv(data, sep="\t", encoding="utf-8")

    refiner = ClusterRefiner(
        threshold = thresh
    )

    refined_df = refiner.refine(
        method = "antismash",
        target_col = "AS_pred"
    )

    print(refined_df)

    refined_df.to_csv(out_file + ".refined.tsv", sep="\t", index=False, header=True)
