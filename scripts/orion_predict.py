import sys
import os
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
ORION = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + "/..")
sys.path.append(ORION)
import pickle
import pandas as pd
import numpy as np
from orion.crf import ClusterCRF
from orion.interface import scripts_interface
from orion.preprocessing import compute_features
from orion.utils import coerce_numeric

SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])))

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_predict.py /home/fleck/scripts/clust/test/test.embed.tsv -o /home/fleck/scripts/clust/test/test

# MAIN
if __name__ == "__main__":
    args = scripts_interface()

    data = args.DATA
    model_file = args.model
    if not model_file:
        model_file = os.path.join(SCRIPT_DIR, "../data/model/feat_v8_param_v2.crf.model")

    out_file = args.out

    e_filter = args.e_filter
    split_col = args.split_col
    sort_cols = args.sort_cols
    weight_col = [coerce_numeric(w) for w in args.w]

    print(args)

    data_tbl = pd.read_csv(data, sep="\t", encoding="utf-8")
    try:
        data_tbl = data_tbl.sort_values(sort_cols).reset_index(drop=True)
    except KeyError as err:
        print("Colums could not be sorted.")
    data_tbl = data_tbl[data_tbl["i_Evalue"] < e_filter]

    for w in weight_col:
        data_tbl = compute_features(data_tbl, weight_type=w)

    print(data_tbl)

    data_tbl.to_csv(out_file + ".format.tsv", sep="\t", index=False, header=True)

    if split_col:
        data_tbl = [s for _, s in data_tbl.groupby(split_col)]
    else:
        data_tbl = [data_tbl]

    with open(model_file, "rb") as f:
        crf = pickle.load(f)

    pred_df = crf.predict_marginals(data=data_tbl)

    pred_df.to_csv(out_file + ".pred.tsv", sep="\t", index=False, header=True)
