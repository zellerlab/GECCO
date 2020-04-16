import typing
from typing import List, Optional

import pandas
from gecco.bgc import Protein, BGC

class ClusterRefiner(object):

    def __init__(
            self,
            threshold: float = 0.4,
            biosynthetic_domains: int = 5,
            seq_col: str ="sequence_id",
            prot_col: str ="protein_id",
            p_col: str ="p_pred",
            domain_col: str = "domain",
            weight_col: str = "log_i_Evalue",
            min_domains: int = 1,
            min_proteins: int = 5,
            join_width: int = 1
    ) -> None:
        self.threshold = threshold
        self.n_biodomains = biosynthetic_domains
        self.n_domains = min_domains
        self.n_proteins = min_proteins
        self.n_proteins = min_proteins
        self.join_width = join_width
        self.seq_col = seq_col
        self.prot_col = prot_col
        self.domain_col = domain_col
        self.weight_col = weight_col
        self.p_col = p_col
        self.grouping = [seq_col, prot_col]

    def find_clusters(
            self,
            domains_df: "pandas.DataFrame",
            method: str = "gecco",
            prefix: str = "cluster",
            lower_threshold: Optional[float] = None,
    ) -> typing.List[BGC]:
        if method == "antismash":
            lt = 0.3 if lower_threshold is None else lower_threshold
        elif method == "gecco":
            lt = self.threshold if lower_threshold is None else lower_threshold
        else:
            raise ValueError(f"unexpected method: {method!r}")
        return self._refine(method, domains_df, lower_threshold=lt)

    def _refine(
        self,
        method: str,
        dataframe: "pandas.DataFrame",
        lower_threshold: float,
    ) -> List[BGC]:
        segments = self.extract_segments(dataframe, lower_threshold)
        clusters = (self._extract_cluster(dataframe, seg) for seg in segments)
        return [bgc for bgc in clusters if bgc.is_valid(criterion=method)]

    def _extract_cluster(
        self,
        dataframe: "pandas.DataFrame",
        segment: "pandas.DataFrame",
    ) -> BGC:
        """Takes a DataFrame and a segement and returns a BGC object"""
        cluster_name = segment["cluster_id"].values[0]
        cluster_prots = set(segment[self.prot_col])
        cluster_df = dataframe.loc[dataframe[self.prot_col].isin(cluster_prots)]
        prot_list = []
        for pid, subdf in segment.groupby(self.prot_col, sort=False):
            protein = Protein(
                seq_id = subdf[self.seq_col].values[0],
                start = subdf["start"].min(),
                end = subdf["end"].max(),
                name = pid,
                domains = subdf[self.domain_col].values,
                weights = subdf[self.weight_col].values,
                probability = subdf[self.p_col].mean(),
            )
            prot_list.append(protein)
        return BGC(prot_list, name=cluster_name)

    def extract_segments(
        self,
        df: "pandas.DataFrame",
        lower_threshold: float,
    ) -> List["pandas.DataFrame"]:
        """
        Extracts segments from a data frame which are determined by p_col.
        Segments are named with prefix_[cluster_number].
        """
        cluster_num = 1
        cluster_state = False
        cluster_list = []
        for n in range(len(df)):
            row = df.iloc[n]
            if row[self.p_col] >= lower_threshold:
                # non-cluster -> cluster
                if not cluster_state:
                    cluster_name = f"{row[self.seq_col]}_cluster_{cluster_num}"
                    row = pandas.DataFrame(row).transpose()
                    cluster_start = row["start"]
                    cluster_df = row
                    cluster_state = True
                # cluster -> cluster
                else:
                    cluster_df = cluster_df.append(row)
                    # Check if last row
                    if n == len(df) - 1:
                        cluster_list.append(
                            cluster_df.assign(idx=n, cluster_id=cluster_name)
                        )
            else:
                # cluster -> non-cluster
                if cluster_state:
                    cluster_list.append(
                        cluster_df.assign(idx=n, cluster_id=cluster_name)
                    )
                    cluster_num += 1
                    cluster_state = False
                # non-cluster -> non-cluster
                # pass
        return cluster_list

    def segment(
        self,
        df: "pandas.DataFrame",
        lower_threshold: float,
        prefix: str,
    ) -> "pandas.DataFrame":
        """
        Determines coordinates of segments determined by p_col over
        a lower_thresh.
        """
        cluster_num = 1
        cluster_state = False
        cluster_list = []
        for n in range(len(df)):
            row = df.iloc[n]
            if row[self.p_col] >= lower_threshold:
                # non-cluster -> cluster
                if not cluster_state:
                    cluster_dict = {}
                    cluster_dict[self.seq_col] = row[self.seq_col]
                    cluster_dict["cluster_id"] = prefix + "_" + str(cluster_num)
                    cluster_dict["start"] = min(row.start, row.end)
                    cluster_state = True
                # cluster -> cluster
                # pass
            else:
                # cluster -> non-cluster
                if cluster_state:
                    cluster_dict["end"] = max(row.start, row.end)
                    cluster_list.append(cluster_dict)
                    cluster_num += 1
                    cluster_state = False
                # non-cluster -> non-cluster
                # pass
        return pandas.DataFrame(cluster_list)
