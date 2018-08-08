import sys
import os
ORION = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + "/..")
sys.path.append(ORION)
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE, MDS
from orion.refine import ClusterRefiner
from orion.interface import scripts_interface
from orion.utils import jsd, jsd_pairwise
from orion.knn import ClusterKNN

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_knn.py /home/fleck/scripts/clust/test/test.pred.tsv -o /home/fleck/scripts/clust/test/test --split-col BGC_id --sort_col BGC_id start -t 1

# MAIN
if __name__ == "__main__":
    # args = scripts_interface()
    #
    # data = args.DATA
    # out_file = args.out
    # thresh = args.thresh
    # split_col = args.split_col
    #
    # print(args)

    data = "~/scripts/clust/test/test.pred.tsv"
    train = ORION + "/data/knn/domain_composition.tsv"
    types = ORION + "/data/knn/type_labels.tsv"
    thresh = 0.5
    split_col = "sequence_id"

    data_df = pd.read_csv(data, sep="\t", encoding="utf-8")

    train_df = pd.read_csv(train, sep="\t", encoding="utf-8")
    comp_array = train_df.iloc[:,1:].values
    id_array = train_df["BGC_id"].values
    pfam_array = train_df.columns.values[1:]

    types_df = pd.read_csv(types, sep="\t", encoding="utf-8")
    types_array = types_df["cluster_type"].values
    subtypes_array = types_df["subtype"].values

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

    cluster_comp = np.array(
        [c.domain_composition(all_possible=pfam_array) for c in cluster_list]
    )

    # pcoa_res = MDS(jsd_mat)
    # tsne = TSNE(metric=jsd_pairwise, perplexity=2, n_iter=2000)
    # tsne_res = tsne.fit_transform(jsd_mat)

    knn = ClusterKNN(metric="jsd", n_neighbors=1)
    knn_pred = knn.fit_predict(comp_array, cluster_comp, y=types_array)

    print(knn_pred)
