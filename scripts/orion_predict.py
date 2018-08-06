import sys
import os
ORION = os.path.abspath(os.path.dirname(sys.argv[0]) + "/..")
sys.path.append(ORION)
import pickle
import pandas as pd
import numpy as np
from itertools import product
from orion.crf import ClusterCRF
from orion.interface import crf_interface
from orion.preprocessing import compute_features
from orion.utils import coerce_numeric

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_predict.py /home/fleck/scripts/clust/test/test.embed.tsv --model /home/fleck/bin/orion/data/model/f5_eval_s_t50.crf.model  -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = crf_interface()

    data = args.DATA
    model_file = args.model
    out_file = args.out

    e_filter = args.e_filter
    split_col = args.split_col
    sort_cols = args.sort_cols
    weight_col = [coerce_numeric(w) for w in args.w]

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    data_tbl = data_tbl.sort_values(sort_cols).reset_index()
    data_tbl = data_tbl[data_tbl["i_Evalue"] < e_filter]

    for w in weight_col:
        data_tbl = compute_features(data_tbl, weight_type=w)

    print(data_tbl)

    data_tbl.to_csv(out_file + ".format.tsv", sep="\t", index=False, header=True)

    data_tbl = [s for _, s in data_tbl.groupby(split_col)]

    with open(model_file, "rb") as f:
        crf = pickle.load(f)

    crf.weights = [1]

    pred_df = crf.predict_marginals(data=data_tbl)

    pred_df.to_csv(out_file + ".pred.tsv", sep="\t", index=False, header=True)
