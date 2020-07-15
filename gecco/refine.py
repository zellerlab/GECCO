"""Algorithm to smooth contiguous BGC predictions into single regions.
"""

import copy
import functools
import operator
import typing
from typing import List, Mapping, Optional, Tuple, Iterator

import pandas
from Bio.SeqRecord import SeqRecord

from .model import Cluster, Domain, Gene, Protein, Strand


# fmt: off
# `set` of `str`: A set of domains from Pfam considered 'biosynthetic' by AntiSMASH.
BIO_PFAMS = frozenset({
    "PF00109", "PF02801", "PF08659", "PF00378", "PF08541", "PF08545",
    "PF02803", "PF00108", "PF02706", "PF03364", "PF08990", "PF00501",
    "PF00668", "PF08415", "PF00975", "PF03061", "PF00432", "PF00494",
    "PF03936", "PF01397", "PF00432", "PF04275", "PF00348", "PF02401",
    "PF04551", "PF00368", "PF00534", "PF00535", "PF02922", "PF01041",
    "PF00128", "PF00908", "PF02719", "PF04321", "PF01943", "PF02806",
    "PF02350", "PF02397", "PF04932", "PF01075", "PF00953", "PF01050",
    "PF03033", "PF01501", "PF05159", "PF04101", "PF02563", "PF08437",
    "PF02585", "PF01721", "PF02052", "PF02674", "PF03515", "PF04369",
    "PF08109", "PF08129", "PF09221", "PF09683", "PF10439", "PF11420",
    "PF11632", "PF11758", "PF12173", "PF04738", "PF04737", "PF04604",
    "PF05147", "PF08109", "PF08129", "PF08130", "PF00155", "PF00202",
    "PF00702", "PF06339", "PF04183", "PF10331", "PF03756", "PF00106",
    "PF01370", "PF00107", "PF08240", "PF00441", "PF02770", "PF02771",
    "PF08028", "PF01408", "PF02894", "PF00984", "PF00725", "PF03720",
    "PF03721", "PF07993", "PF02737", "PF00903", "PF00037", "PF04055",
    "PF00171", "PF00067", "PF01266", "PF01118", "PF02668", "PF00248",
    "PF01494", "PF01593", "PF03992", "PF00355", "PF01243", "PF00384",
    "PF01488", "PF00857", "PF04879", "PF08241", "PF08242", "PF00698",
    "PF00483", "PF00561", "PF00583", "PF01636", "PF01039", "PF00288",
    "PF00289", "PF02786", "PF01757", "PF02785", "PF02409", "PF01553",
    "PF02348", "PF00891", "PF01596", "PF04820", "PF02522", "PF08484",
    "PF08421",
})


class ClusterRefiner:
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
    ) -> Iterator[Cluster]:
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
        _validate_cluster = functools.partial(self._validate_cluster, criterion=criterion)
        lt = self.threshold if lower_threshold is None else lower_threshold
        clusters = map(_extract_cluster, self._iter_segments(features, lt))
        return filter(_validate_cluster, clusters)

    def _extract_cluster(
        self,
        genes: Mapping[str, Gene],
        segment: "pandas.DataFrame",
    ) -> Cluster:
        """Take a segment of contiguous domains and returns a `BGC`.
        """
        cluster_genes = []
        for prot_id, subdf in segment.groupby(self.protein_column, sort=False):
            gene = copy.deepcopy(genes[prot_id])
            gene.probability = subdf[self.probability_column].mean()
            gene.protein.domains = [
                Domain(t.domain, t.domain_start, t.domain_end, t.i_Evalue)
                for t in subdf.itertuples()
            ]
            cluster_genes.append(gene)
        return Cluster(
            id=segment.cluster_id.values[0],
            genes=cluster_genes,
        )

    def _validate_cluster(
        self,
        cluster: Cluster,
        criterion: str = "gecco",
        threshold: float = 0.6,
        n_cds: int = 5,
    ) -> bool:
        if criterion == "gecco":
            return True
        elif criterion == "antismash":
            domains = {d.name for gene in cluster.genes for d in gene.protein.domains}
            p_crit = numpy.mean([g.probability for g in cluster.genes]) >= threshold
            bio_crit = len(domains & BIO_PFAMS) >= n_biopfams
            cds_crit = len(cluster.genes) >= n_cds
            return p_crit # TODO: & bio_crit & cds_crit
        else:
            raise ValueError(f"unknown BGC filtering criterion: {criterion}")

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
