"""Algorithm to smooth contiguous BGC predictions into single regions.
"""

import functools
import operator
import typing
from typing import List, Mapping, Optional, Tuple, Iterator

import pandas
from Bio.SeqRecord import SeqRecord

from .bgc import Protein, BGC
from .model import Cluster, ClusterGene, ClusterProtein, Domain, Gene, Strand


class ClusterRefiner(object):
    """A post-processor to extract contiguous BGCs from CRF predictions.

    Like `~gecco.crf.ClusterCRF`, it can be configured to support non-standard
    column names in arguments `~pandas.DataFrame`.

    """

    def __init__(
        self,
        threshold: float = 0.4,
        sequence_column: str = "sequence_id",
        protein_column: str = "protein_id",
        probability_column: str = "p_pred",
        domain_column: str = "domain",
        weight_column: str = "rev_i_Evalue",
    ) -> None:
        """Create a new `ClusterRefiner` instance.

        Arguments:
            threshold (`float`): The probability threshold to use to consider
                a protein to be part of a BGC region.
            sequence_column (`str`): The name of the column containing the
                sequence ID for each row.
            protein_column (`str`): The name of the column containing the
                protein ID for each row.
            probability_column (`str`): The name of the columnn containing BGC
                probability predicted by `~gecco.crf.ClusterCRF`.
            domain_column (`str`): The name of the column containing the domain
                name for each row.
            weight_column (`str`): The name of the column containing the weight
                for each row.

        """
        self.threshold = threshold
        self.sequence_column = sequence_column
        self.protein_column = protein_column
        self.probability_column = probability_column
        self.domain_column = domain_column
        self.weight_column = weight_column

    def iter_clusters(
        self,
        genes: Mapping[str, Gene],
        features: "pandas.DataFrame",
        criterion: str = "gecco",
        prefix: str = "cluster",
        lower_threshold: Optional[float] = None,
    ) -> Iterator[BGC]:
        """Find all clusters in a table of CRF predictions.

        Arguments:
            features (`~pandas.DataFrame`): The data frame containing BGC
                probability predictions, obtained by `~gecco.crf.ClusterCRF`.
            criterion (`str`): The criterion to use when checking for BGC
                validity. See `gecco.bgc.BGC.is_valid` documentation for
                allowed values and expected behaviours.
            prefix (`str`): The name prefix to use to label BGCs. Each BGC is
                labelled ``{prefix}_{index}``, with ``index`` starting at 1.
            lower_threshold (`float`, optional): If given, overrides the
                object probability threshold (``self.threshold``) to extract
                BGC regions.

        Yields:
            `BGC`: Contiguous biosynthetic gene clusters found the input.

        Raises:
            `ValueError`: When ``criterion`` is not an allowed value.

        """
        _extract_cluster = functools.partial(self._extract_cluster, genes)
        lt = self.threshold if lower_threshold is None else lower_threshold
        clusters = map(_extract_cluster, self._iter_segments(features, lt))
        return filter(lambda bgc: bgc.is_valid(criterion=criterion), clusters)

    def _extract_cluster(
        self,
        genes: Mapping[str, Gene],
        segment: "pandas.DataFrame",
    ) -> BGC:
        """Take a segment of contiguous domains and returns a `BGC`.
        """
        cluster_genes = []
        for prot_id, subdf in segment.groupby(self.protein_column, sort=False):
            domains = [
                Domain(row.domain, row.domain_start, row.domain_end, row.i_Evalue)
                for row in subdf.itertuples()
            ]
            cluster_protein = ClusterProtein(
                id=prot_id,
                seq=genes[prot_id].protein.seq,
                domains=domains,
            )
            cluster_gene = ClusterGene(
                seq_id=genes[prot_id].seq_id,
                start=genes[prot_id].start,
                end=genes[prot_id].end,
                strand=genes[prot_id].strand,
                protein=cluster_protein,
            )
            cluster_genes.append(cluster_gene)
        return Cluster(
            id=segment.cluster_id.values[0],
            genes=cluster_genes,
        )

        # return BGC(
        #     name=segment.cluster_id.values[0],
        #     proteins=[
        #         Protein(
        #             seq_id=prot_df[self.sequence_column].values[0],
        #             start=int(prot_df.start.min()),
        #             end=int(prot_df.end.max()),
        #             name=prot_id,
        #             strand=prot_df.strand.values[0],
        #             domains=prot_df[self.domain_column].values,
        #             weights=prot_df[self.weight_column].values,
        #             probability=prot_df[self.probability_column].mean(),
        #         )
        #         for prot_id, prot_df in segment.groupby(self.protein_column, sort=False)
        #     ],
        # )

    def _iter_segments(
        self,
        features: "pandas.DataFrame",
        lower_threshold: float,
    ) -> Iterator["pandas.DataFrame"]:
        """Iterate over contiguous BGC segments from a CRF prediction table.

        Yields:
            `~pandas.DataFrame`: A data frame containing one or more
            consecutive proteins, found to be part of the same gene cluster.

        """
        in_cluster: int = False
        cluster_num: int = 1
        cluster: List[Tuple[...]]

        # getter functions to get the probability and the seq_id from the
        # namedtuple fields, which do not support item getter protocol
        p_getter = operator.attrgetter(self.probability_column)
        seq_getter = operator.attrgetter(self.sequence_column)

        # process all proteins to find potential BGCs
        for idx, row in enumerate(features.itertuples()):
            if p_getter(row) >= lower_threshold:
                # non-cluster -> cluster
                if not in_cluster:
                    id_ = f"{seq_getter(row)}_cluster_{cluster_num}"
                    cluster = [row]
                    in_cluster = True
                # cluster -> cluster
                else:
                    cluster.append(row)
            elif in_cluster:
                # cluster -> non-cluster
                yield pandas.DataFrame(cluster).assign(idx=idx, cluster_id=id_)
                cluster_num += 1
                in_cluster = False
                # non-cluster -> non-cluster
                # do nothing #

        # if the sequence is over while we were still in a cluster, we need
        # to recover the last cluster
        if in_cluster:
            yield pandas.DataFrame(cluster).assign(idx=idx, cluster_id=id_)
