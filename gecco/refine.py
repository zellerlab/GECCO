"""Algorithm to smooth contiguous gene cluster predictions into single regions.
"""

import collections.abc
import itertools
import functools
import operator
import typing
from typing import List, Mapping, Optional, Tuple, Iterator

import numpy
from Bio.SeqRecord import SeqRecord

from .model import Cluster, Domain, Gene, Protein, Strand


__all__ = ["BIO_PFAMS", "GeneGrouper", "ClusterRefiner"]


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



class GeneGrouper:
    """A callable to group genes under or over a probability threshold.

    Use with a list of genes in combination with `itertools.groupby`.
    """

    def __init__(self, threshold: float):  # noqa: D102, D107
        self.in_cluster = False
        self.threshold = threshold

    def __call__(self, gene: Gene) -> bool:  # noqa: D102
        if gene.average_probability is not None:
            self.in_cluster = gene.average_probability > self.threshold
        return self.in_cluster



class ClusterRefiner:
    """A post-processor to extract contiguous clusters from CRF predictions.
    """

    def __init__(
        self,
        threshold: float = 0.8,
        criterion: str = "gecco",
        n_cds: int = 5,
        n_biopfams: int = 5,
        average_threshold: float = 0.6,
        edge_distance: int = 0,
    ) -> None:
        """Create a new `ClusterRefiner` instance.

        Arguments:
            threshold (`float`): The probability threshold to use to consider
                a protein to be part of a gene cluster.
            criterion (`str`): The criterion to use when checking for cluster
                validity. 
            n_cds (`int`): The minimum number of genes a gene cluster must
                contain to be considered valid. If ``criterion`` is ``gecco``,
                then this is the minimum number of **annotated** CDS.
            n_biopfams (`int`): The minimum number of biosynthetic Pfam
                domains a gene cluster must contain to be considered valid
                (*only when the criterion is* ``antismash``).
            average_threshold (`int`): The average probability threshold to
                use to consider a gene cluster valid (*only when the 
                criterion is* ``antismash``).
            edge_distance (`int`): The minimum distance from the edge the 
                gene cluster must be located (it may start at an edge, but must 
                span for longer than ``edge_distance``), in number of annotated 
                genes (*only when the criterion is* ``gecco``).

        """
        self.threshold = threshold
        self.criterion = criterion
        self.n_cds = n_cds
        self.n_biopfams = n_biopfams
        self.average_threshold = average_threshold
        self.edge_distance = edge_distance

    def iter_clusters(self, genes: List[Gene]) -> Iterator[Cluster]:
        """Find all clusters in a table of CRF predictions.

        Arguments:
            genes (`list` of `~gecco.model.Gene`): A list of genes with
                probability annotations estimated by `~gecco.crf.ClusterCRF`.

        Yields:
            `gecco.model.Cluster`: Valid clusters found in the input with
            respect to the postprocessing criterion given at initialisation.

        """
        for seq, cluster in self._iter_clusters(genes):
            trimmed = self._trim_cluster(cluster)
            if self._validate_cluster(seq, cluster):
                yield cluster
        #
        # unfiltered_clusters = map(self._trim_cluster, self._iter_clusters(genes))
        # return filter(self._validate_cluster, unfiltered_clusters)

    def _validate_cluster(self, seq: List[Gene], cluster: Cluster) -> bool:
        """Check a cluster validity depending on the postprocessing criterion.
        """
        if self.criterion == "gecco":
            # check that the number of cluster genes that are outside of
            # the edge region is above the number of required CDS
            annotated = [g for g in cluster.genes if g.protein.domains]
            cds_crit = len(annotated) >= self.n_cds
            # extract IDs of annotated genes that are too close to the edge
            if self.edge_distance > 0:
                annotated_ids = [g.id for g in seq if g.protein.domains]
                edge_genes = set(annotated_ids[:self.edge_distance]).union(annotated_ids[-self.edge_distance:])
            else:
                edge_genes = set()
            # NOTE (@althonos): it is needed for compatibility with the post-processed results
            #                   that we filter on any number of cluster genes, but it would be
            #                   better to filter on the number of annotated cluster genes instead.
            edge_crit = len(set(g.id for g in cluster.genes).difference(edge_genes)) >= self.n_cds
            return cds_crit and edge_crit
        elif self.criterion == "antismash":
            domains = {d.name for gene in cluster.genes for d in gene.protein.domains}
            p_crit = numpy.mean([typing.cast(float, g.average_probability) for g in cluster.genes]) >= self.average_threshold
            bio_crit = len(domains & BIO_PFAMS) >= self.n_biopfams
            cds_crit = len(cluster.genes) >= self.n_cds
            return p_crit and bio_crit and cds_crit
        else:
            raise ValueError(f"Unknown cluster filtering criterion: {self.criterion}")

    def _trim_cluster(self, cluster: Cluster) -> Cluster:
        """Remove unannotated proteins from the cluster edges.
        """
        while cluster.genes and not cluster.genes[0].protein.domains:
            cluster.genes.pop(0)
        while cluster.genes and not cluster.genes[-1].protein.domains:
            cluster.genes.pop()
        return cluster

    def _iter_clusters(
        self,
        genes: List[Gene],
    ) -> Iterator[Tuple[List[Gene], Cluster]]:
        """Iterate over contiguous cluster segments from a list of genes.
        """
        grouper = GeneGrouper(self.threshold)
        key = operator.attrgetter("source.id")

        # iterate over the genes grouped by sequence ids
        for seq_id, sequence in itertools.groupby(sorted(genes, key=key), key=key):
            # sort genes within the same sequence by coordinates
            seqsort = sorted(sequence, key=operator.attrgetter("start", "end"))
            # group contiguous genes if they are over or under the threshold
            groups = itertools.groupby(seqsort, key=grouper)
            # filter out regions that are not identified to be clusters
            clusters = (genes for in_cluster, genes in groups if in_cluster)
            for i, cluster in enumerate(clusters):
                yield seqsort, Cluster(id=f"{seq_id}_cluster_{i+1}", genes=list(cluster))
