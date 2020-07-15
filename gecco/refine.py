"""Algorithm to smooth contiguous BGC predictions into single regions.
"""

import itertools
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
        criterion: str = "gecco",
        n_cds: int = 5,
        n_biopfams: int = 5,
        average_threshold: float = 0.6,
    ) -> None:
        """Create a new `ClusterRefiner` instance.

        Arguments:
            threshold (`float`): The probability threshold to use to consider
                a protein to be part of a BGC region.
            criterion (`str`): The criterion to use when checking for BGC
                validity. See `gecco.bgc.BGC.is_valid` documentation for
                allowed values and expected behaviours.

        """
        self.threshold = threshold
        self.criterion = criterion
        self.n_cds = n_cds
        self.n_biopfams = n_biopfams
        self.average_threshold = average_threshold

    def iter_clusters(self, genes: List[Gene]) -> Iterator[Cluster]:
        """Find all clusters in a table of CRF predictions.

        Arguments:
            genes (`list` of `~gecco.model.Gene`): A list of genes with
                probability annotations estimated by `~gecco.crf.ClusterCRF`.

        Yields:
            `gecco.model.Cluster`: Valid clusters found in the input with
            respect to the postprocessing criterion given at initialisation.

        """
        return filter(self._validate_cluster, self._iter_clusters(genes))

    def _validate_cluster(self, cluster: Cluster) -> bool:
        if self.criterion == "gecco":
            cds_crit = len(cluster.genes) >= self.n_cds
            return cds_crit
        elif self.criterion == "antismash":
            domains = {d.name for gene in cluster.genes for d in gene.protein.domains}
            p_crit = numpy.mean([g.probability for g in cluster.genes]) >= self.average_threshold
            bio_crit = len(domains & BIO_PFAMS) >= self.n_biopfams
            cds_crit = len(cluster.genes) >= self.n_cds
            return p_crit & bio_crit & cds_crit
        else:
            raise ValueError(f"unknown BGC filtering criterion: {self.criterion}")

    def _iter_clusters(
        self,
        genes: List[Gene],
    ) -> Iterator[Cluster]:
        """Iterate over contiguous BGC segments from a list of genes.
        """
        thr = self.threshold
        k1 = operator.attrgetter("seq_id")
        k2 = operator.attrgetter("start", "end")

        for seq_id, sequence in itertools.groupby(sorted(genes, key=k1), key=k1):
            cluster_num = 1
            in_cluster = False

            for gene in sorted(sequence, key=k2):
                if not in_cluster:
                    # not cluster -> cluster
                    if gene.probability is not None and gene.probability >= thr:
                        cluster = Cluster(id=f"{gene.seq_id}_cluster_{cluster_num}")
                        cluster.genes = [gene]
                        in_cluster = True
                    # not cluster -> not cluster
                    # pass
                else:
                    # cluster -> cluster
                    if gene.probability is None or gene.probability >= thr:
                        cluster.genes.append(gene)
                    # cluster -> not cluster
                    else:
                        yield cluster
                        in_cluster = False
                        cluster_num += 1

            # if the sequence is over while we were still in a cluster, we need
            # to recover the last cluster
            if in_cluster:
                yield cluster
