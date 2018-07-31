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

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_predict.py /home/fleck/scripts/clust/test/test.embed.tsv --model /home/fleck/bin/orion/data/model/f5_eval_s_t50.crf.model  -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = crf_interface()

    data = args.DATA
    model_file = args.model
    out_file = args.out

    C1 = 0.15
    C2 = 1.75

    weight_col = args.w
    feature_type = args.feature_type
    e_filter = args.e_filter
    split_col = args.split_col
    overlap = args.overlap

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    data_tbl = data_tbl.sort_values(["genome_id", "start", "domain_start"])

    # data_tbl = data_tbl[data_tbl["i_Evalue"] < e_filter]
    # data_tbl = compute_features(data_tbl, split_col=split_col)
    #
    print(data_tbl)
    #
    # data_tbl.to_csv(out_file + ".format.tsv", sep="\t", index=False, header=True)

    data_tbl = [s for _, s in data_tbl.groupby(split_col)]


    crf = ClusterCRF(data_tbl,
        feature_cols = ["pfam"],
        group_col = "protein_id",
        weight_cols = [weight_col],
        feature_type = feature_type,
        overlap = overlap,
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2
    )

    with open(model_file, "rb") as f:
        model = pickle.load(f)
        crf.model = model

    pred_df = crf.predict_marginals()

    print(pred_df)

    pred_df.to_csv(out_file + ".pred.tsv", sep="\t", index=False, header=True)