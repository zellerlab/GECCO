import sys
import os
ORION = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + "/..")
sys.path.append(ORION)
import numpy as np
import pandas as pd
from skbio.stats.ordination import pcoa
from sklearn.manifold import TSNE
from sklearn.neighbors import KNeighborsClassifier
from orion.refine import ClusterRefiner
from orion.interface import scripts_interface
from orion.utils import jsd

### TEST ###
# python /home/fleck/bin/orion/scripts/orion_knn.py /home/fleck/scripts/clust/test/test.pred.tsv -o /home/fleck/scripts/clust/test/test

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

    print(cluster_comp)

    jsd_mat = jsd(cluster_comp)

    print(jsd_mat)

    pcoa_res = pcoa(jsd_mat)
    tsne = TSNE(metric="precomputed", perplexity=2, n_iter=2000)
    tsne_res = tsne.fit_transform(jsd_mat)

    knn = KNeighborsClassifier(n_neighbors=1, metric="precomputed")
    knn.fit(jsd_mat, y=["bla", "blubb", "bla"])
    knn_pred = knn.predict(jsd_mat)

    print(knn_pred)
