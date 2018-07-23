import random
import argparse
import pandas as pd
import numpy as np
import multiprocessing
from itertools import product
from crf import ClusterCRF

### TEST ###
# python /home/fleck/scripts/clust/main/crf_cvpy /home/fleck/scripts/clust/test/test.embed.tsv -o /home/fleck/scripts/clust/test/test


# FUNC
def interface():
    parser = argparse.ArgumentParser(description="Takes a set of training and cross- validation samples and writes the CV results as tsv to file.")

    parser.add_argument("DATA",
                        type=str,
                        metavar="<DATA>",
                        help="Pfam tables with training instances.")

    parser.add_argument("-o", "--output-basename",
                        dest="out",
                        type=str,
                        default="CRF",
                        metavar="<basename>",
                        help="Basename for the output table with the CV results.")

    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        metavar="<threads>",
                        help="Number of CPUs.")

    parser.add_argument("-e", "--evalue",
                        dest="e_filter",
                        type=float,
                        default="1",
                        help="E-value threshold for the test set.")

    parser.add_argument("-w", "--weight-col",
                        dest="w",
                        default="pseudo_norm_prot",
                        type=str,
                        help="Column to be used as local weights on pfam domains.")

    parser.add_argument("--feature-type",
                        dest="feature_type",
                        type=str,
                        default="group",
                        help="How features should be extracted. 'Single', 'overlap' or on some grouping level ('group').")

    parser.add_argument("--truncate",
                        dest="truncate",
                        type=int,
                        help="Training set will be truncated to this length.")

    parser.add_argument("--overlap",
                        dest="overlap",
                        type=int,
                        default="2",
                        help="If overlapping features: How much overlap.")

    parser.add_argument("--no-shuffle",
                        dest="shuffle",
                        action="store_false",
                        help="Switch to turn of shuffling of the data before doing CV.")

    args = parser.parse_args()
    return args


# MAIN
if __name__ == "__main__":
    args = interface()

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

    crf = ClusterCRF(data_tbl, "BGC",
        feature_cols = ["pfam"],
        weight_cols = [weight_col],
        feature_type = feature_type,
        overlap = overlap,
        algorithm = "lbfgs",
        c1 = C1,
        c2 = C2)

    results = crf.cv(
        threads = threads,
        e_filter = e_filter,
        truncate = truncate,
        strat_col = "BGC_type")

    result_df = (pd .concat(results)
                    .assign(c1 = C1,
                        c2 = C2,
                        feature_type = feature_type,
                        e_filter = e_filter,
                        overlap = overlap,
                        weight = weight_col,
                        truncate = truncate,
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
