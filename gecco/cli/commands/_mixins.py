import collections
import csv
import gzip
import io
import itertools
import math
import multiprocessing.pool
import os
import operator
import typing
from typing import (
    BinaryIO,
    Container,
    Collection,
    Iterator,
    Iterable,
    Optional,
    List,
    Union
)

from ..._meta import zopen
from .._utils import ProgressReader, guess_sequences_format
from ._base import Command, CommandExit, InvalidArgument

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord
    from ...model import Cluster, Gene, FeatureTable, GeneTable, ClusterTable
    from ...hmmer import HMM
    from ...types import TypeClassifier


class SequenceLoaderMixin(Command):

    format: Optional[str]
    genome: str

    def _load_sequences(self) -> Iterator["SeqRecord"]:
        from Bio import SeqIO

        try:
            # guess format or use the one given in CLI
            if self.format is not None:
                format: Optional[str] = self.format.lower()
                self.info("Using", "user-provided sequence format", repr(format), level=2)
            else:
                self.info("Detecting", "sequence format from file contents", level=2)
                format = guess_sequences_format(self.genome)
                if format is None:
                    raise RuntimeError(f"Failed to detect format of {self.genome!r}")
                self.success("Detected", "format of input as", repr(format), level=2)
            # get filesize and unit
            input_size = os.stat(self.genome).st_size
            total, scale, unit = ProgressReader.scale_size(input_size)
            task = self.progress.add_task("Loading sequences", total=total, unit=unit, precision=".1f")
            # load sequences
            n = 0
            self.info("Loading", "sequences from genomic file", repr(self.genome), level=1)
            with zopen(self.genome) as file:
                with ProgressReader(file, self.progress, task, scale) as reader:
                    for record in SeqIO.parse(io.TextIOWrapper(reader), format):  # type: ignore
                        yield record
                        n += 1
        except FileNotFoundError as err:
            self.error("Could not find input file:", repr(self.genome))
            raise CommandExit(err.errno) from err
        except ValueError as err:
            self.error("Failed to load sequences:", err)
            raise CommandExit(getattr(err, "errno", 1)) from err
        else:
            self.success("Found", n, "sequences", level=1)


class TableLoaderMixin(Command):

    genes: str
    features: List[str]

    def _load_genes(self) -> Iterator["Gene"]:
        from ...model import GeneTable

        try:
            # get filesize and unit
            input_size = os.stat(self.genes).st_size
            total, scale, unit = ProgressReader.scale_size(input_size)
            task = self.progress.add_task("Loading genetable", total=total, unit=unit, precision=".1f")
            # load gene table
            self.info("Loading", "genes table from file", repr(self.genes))
            with typing.cast(BinaryIO, ProgressReader(open(self.genes, "rb"), self.progress, task, scale)) as in_:
                table = GeneTable.load(in_)
            # count genes and yield gene objects
            n_genes = len(set(table.protein_id))
            unit = "gene" if n_genes == 1 else "genes"
            task = self.progress.add_task("Building genes", total=n_genes, unit=unit, precision="")
            yield from self.progress.track(table.to_genes(), task_id=task)
        except OSError as err:
            self.error("Fail to parse genes coordinates: {}", err)
            raise CommandExit(err.errno) from err

    def _load_features(self) -> "FeatureTable":
        from ...model import FeatureTable

        features = FeatureTable()
        for filename in self.features:
            try:
                # get filesize and unit
                input_size = os.stat(filename).st_size
                total, scale, unit = ProgressReader.scale_size(input_size)
                task = self.progress.add_task("Loading features", total=total, unit=unit, precision=".1f")
                # load features
                self.info("Loading", "features table from file", repr(filename))
                with zopen(filename) as file:
                    with ProgressReader(file, self.progress, task, scale) as in_:
                        features += FeatureTable.load(in_)
            except FileNotFoundError as err:
                self.error("Could not find feature file:", repr(filename))
                raise CommandExit(err.errno) from err

        self.success("Loaded", "a total of", len(features), "features", level=1)
        return features

    def _annotate_genes(self, genes: List["Gene"], features: "FeatureTable") -> List["Gene"]:
        from ...model import Domain

        # index genes by protein_id
        gene_index = { gene.protein.id:gene for gene in genes }
        if len(gene_index) < len(genes):
            raise ValueError("Duplicate gene names in input genes")

        # add domains from the feature table
        unit = "row" if len(features) == 1 else "rows"
        task = self.progress.add_task("Annotating genes", total=len(features), unit=unit, precision="")
        for i in self.progress.track(range(len(features)), task_id=task):
            # get gene by ID and check consistency
            length = features.end[i] - features.start[i]
            gene = gene_index[features.protein_id[i]]
            if gene.source.id != features.sequence_id[i]:
                raise ValueError(f"Mismatched source sequence for {features.protein_id[i]!r}: {gene.source.id!r} != {features.sequence_id[i]!r}")
            elif gene.end - gene.start != length:
                raise ValueError(f"Mismatched gene length for {row.protein_id!r}: {gene.end - gene.start!r} != {length!r}")
            elif gene.start != features.start[i]:
                raise ValueError(f"Mismatched gene start for {row.protein_id!r}: {gene.start!r} != {features.start[i]!r}")
            elif gene.end != features.end[i]:
                raise ValueError(f"Mismatched gene end for {row.protein_id!r}: {gene.end!r} != {features.end[i]!r}")
            elif gene.strand.sign != features.strand[i]:
                raise ValueError(f"Mismatched gene strand for {row.protein_id!r}: {gene.strand.sign!r} != {features.strand[i]!r}")
            # add the row domain to the gene
            domain = Domain(
                name=features.domain[i],
                start=features.domain_start[i],
                end=features.domain_end[i],
                hmm=features.hmm[i],
                i_evalue=features.i_evalue[i],
                pvalue=features.pvalue[i],
                probability=features.cluster_probability[i],
            )
            gene.protein.domains.append(domain)

        # return
        return list(gene_index.values())


class OutputWriterMixin(Command):

    genome: str
    output_dir: str

    def _make_output_directory(self, outputs: List[str]) -> None:
        # Make output directory
        self.info("Using", "output folder", repr(self.output_dir), level=1)
        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except OSError as err:
            self.error("Could not create output directory: {}", err)
            raise CommandExit(err.errno) from err

        # Check if output files already exist
        for output in outputs:
            if os.path.isfile(os.path.join(self.output_dir, output)):
                self.warn("Output folder contains files that will be overwritten")
                break

    def _write_feature_table(self, genes: List["Gene"]) -> None:
        from ...model import FeatureTable

        base, _ = os.path.splitext(os.path.basename(self.genome))
        pred_out = os.path.join(self.output_dir, f"{base}.features.tsv")
        self.info("Writing", "feature table to", repr(pred_out), level=1)
        with open(pred_out, "wb") as f:
            FeatureTable.from_genes(genes).dump(f)

    def _write_genes_table(self, genes: List["Gene"]) -> None:
        from ...model import GeneTable

        base, _ = os.path.splitext(os.path.basename(self.genome))
        pred_out = os.path.join(self.output_dir, f"{base}.genes.tsv")
        self.info("Writing", "gene table to", repr(pred_out), level=1)
        with open(pred_out, "wb") as f:
            GeneTable.from_genes(genes).dump(f)

    def _write_cluster_table(self, clusters: List["Cluster"]) -> None:
        from ...model import ClusterTable

        base, _ = os.path.splitext(os.path.basename(self.genome))
        cluster_out = os.path.join(self.output_dir, f"{base}.clusters.tsv")
        self.info("Writing", "cluster table to", repr(cluster_out), level=1)
        with open(cluster_out, "wb") as out:
            ClusterTable.from_clusters(clusters).dump(out)

    def _write_clusters(self, clusters: List["Cluster"], merge: bool = False) -> None:
        from Bio import SeqIO

        if merge:
            base, _ = os.path.splitext(os.path.basename(self.genome))
            gbk_out = os.path.join(self.output_dir, f"{base}.clusters.gbk")
            records = (cluster.to_seq_record() for cluster in clusters)
            SeqIO.write(records, gbk_out, "genbank")
        else:
            for cluster in clusters:
                gbk_out = os.path.join(self.output_dir, f"{cluster.id}.gbk")
                self.info("Writing", f"cluster [bold blue]{cluster.id}[/] to", repr(gbk_out), level=1)
                SeqIO.write(cluster.to_seq_record(), gbk_out, "genbank")


class DomainFilterMixin(Command):

    e_filter: Optional[float]
    p_filter: Optional[float]
    disentangle: bool

    def _filter_domains(self, genes: List["Gene"]) -> List["Gene"]:
        # Filter i-evalue and p-value if required
        if self.e_filter is not None:
            self.info("Excluding", "domains with e-value over", self.e_filter, level=1)
            key = lambda d: d.i_evalue < self.e_filter
            genes = [
                gene.with_protein(gene.protein.with_domains(filter(key, gene.protein.domains)))
                for gene in genes
            ]
        if self.p_filter is not None:
            self.info("Excluding", "domains with p-value over", self.p_filter, level=1)
            key = lambda d: d.pvalue < self.p_filter
            genes = [
                gene.with_protein(gene.protein.with_domains(filter(key, gene.protein.domains)))
                for gene in genes
            ]
        if self.p_filter is not None or self.e_filter is not None:
            count = sum(1 for gene in genes for domain in gene.protein.domains)
            self.info("Using", "remaining", count, "domains", level=1)
        return genes

    def _disentangle(self, gene: "Gene") -> "Gene":
        if len(gene.protein.domains) <= 1:
            return gene

        domains_to_keep = []
        gene_domains = gene.protein.domains.copy()
        while gene_domains:
            # get a domain and collect other overlaping domains
            domain = gene_domains.pop()
            overlaps = [
                other for other in gene_domains
                if other.start <= domain.end and domain.start <= other.end
            ]
            # keep only overlapping domain with the best p-value
            if not overlaps or domain.pvalue < min(d.pvalue for d in overlaps):
                domains_to_keep.append(domain)
                for other in overlaps:
                    gene_domains.remove(other)

        return gene.with_protein(gene.protein.with_domains(domains_to_keep))

    def _disentangle_domains(self, genes: List["Gene"]) -> List["Gene"]:
        self.info("Disentangling", "overlapping domains in each gene", level=1)
        return list(map(self._disentangle, genes))


class AnnotatorMixin(DomainFilterMixin):

    hmm: Optional[List[str]]
    jobs: int
    bit_cutoffs: Optional[str]

    def _custom_hmms(self) -> Iterable["HMM"]:
        from ...hmmer import HMM

        for path in typing.cast(List[str], self.hmm):
            base = os.path.basename(path)
            file: BinaryIO = zopen(path)
            if base.endswith((".gz", ".lz4", ".xz", ".bz2")):
                base, _ = os.path.splitext(base)
            base, _ = os.path.splitext(base)
            yield HMM(
                id=base,
                version="?",
                url="?",
                path=path,
                size=None,
                relabel_with=r"s/([^\.]*)(\..*)?/\1/"
            )

    def _annotate_domains(self, genes: List["Gene"], whitelist: Optional[Collection[str]] = None) -> List["Gene"]:
        from ...hmmer import PyHMMER, embedded_hmms

        self.info("Running", "HMMER domain annotation", level=1)

        # Run all HMMs over ORFs to annotate with protein domains
        hmms = list(self._custom_hmms() if self.hmm else embedded_hmms())
        total = None if whitelist is None else len(whitelist)
        task_hmms = self.progress.add_task(description=f"Annotating domains", unit="HMMs", total=len(hmms), precision="")
        task_domains = self.progress.add_task(description="", total=total, unit="domains", precision="")
        for hmm in self.progress.track(hmms, task_id=task_hmms, total=len(hmms)):
            self.progress.update(task_domains, description=f" {hmm.id} v{hmm.version}")
            callback = lambda h, t: self.progress.update(task_domains, advance=1)
            self.info("Starting", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
            genes = PyHMMER(hmm, self.jobs, whitelist).run(genes, progress=callback, bit_cutoffs=self.bit_cutoffs)
            self.success("Finished", f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]", level=2)
        self.progress.update(task_id=task_domains, visible=False)

        # Count number of annotated domains
        count = sum(1 for gene in genes for domain in gene.protein.domains)
        self.success("Found", count, "domains across all proteins", level=1)

        # Disentangle if required
        if self.disentangle:
            genes = self._disentangle_domains(genes)

        # Filter i-evalue and p-value if required
        genes = self._filter_domains(genes)

        # Sort genes
        self.info("Sorting", "genes by coordinates", level=2)
        genes.sort(key=lambda g: (g.source.id, g.start, g.end))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

        return genes


class PredictorMixin(Command):
    """Common code for a command that runs probability prediction.
    """

    model: Optional[str]
    no_pad: bool
    threshold: float
    postproc: str
    cds: int
    edge_distance: int

    def _predict_probabilities(self, genes: List["Gene"]) -> List["Gene"]:
        from ...crf import ClusterCRF

        if self.model is None:
            self.info("Loading", "embedded CRF pre-trained model", level=1)
        else:
            self.info("Loading", "CRF pre-trained model from", repr(self.model), level=1)
        model = ClusterCRF.trained(self.model)

        self.info("Predicting", "cluster probabilitites with the model", level=1)
        unit = "batches" if len(genes) > 1 else "batch"
        task = self.progress.add_task("Predicting marginals", total=None, unit=unit, precision="")

        def progress_callback(i: int, total: int) -> None:
            self.progress.update(task_id=task, completed=i, total=total)

        return model.predict_probabilities(
            genes,
            pad=not self.no_pad,
            progress=progress_callback,
        )

    def _extract_clusters(self, genes: List["Gene"]) -> List["Cluster"]:
        from ...refine import ClusterRefiner

        self.info("Extracting", "predicted clusters", level=1)
        refiner = ClusterRefiner(
            threshold=self.threshold,
            criterion=self.postproc,
            n_cds=self.cds,
            edge_distance=self.edge_distance
        )

        total = len({gene.source.id for gene in genes})
        unit = "contigs" if total > 1 else "contig"
        task = self.progress.add_task("Extracting clusters", total=total, unit=unit, precision="")

        clusters: List["Cluster"] = []
        gene_groups = itertools.groupby(genes, lambda g: g.source.id)  # type: ignore
        for _, gene_group in self.progress.track(gene_groups, task_id=task, total=total):
            clusters.extend(refiner.iter_clusters(list(gene_group)))

        return clusters

    def _load_type_classifier(self) -> "TypeClassifier":
        from ...types import TypeClassifier

        self.info("Loading", "type classifier from internal model", level=2)
        return TypeClassifier.trained(self.model)

    def _predict_types(self, clusters: List["Cluster"], classifier: "TypeClassifier") -> List["Cluster"]:
        self.info("Predicting", "gene cluster types", level=1)

        unit = "cluster" if len(clusters) == 1 else "clusters"
        task = self.progress.add_task("Predicting types", total=len(clusters), unit=unit, precision="")

        clusters_new = []
        for cluster in self.progress.track(clusters, task_id=task):
            clusters_new.extend(classifier.predict_types([cluster]))
            if cluster.type:
                name = "/".join(f"[bold blue]{name}[/]" for name in cluster.type.names)
                prob = "/".join(f"[bold purple]{cluster.type_probabilities[name]:.0%}[/]" for name in cluster.type.names)
                self.success(f"Predicted type of [bold blue]{cluster.id}[/] as {name} ({prob} confidence)")
            else:
                ty = max(cluster.type_probabilities, key=cluster.type_probabilities.get)   # type: ignore
                prob = f"[bold purple]{cluster.type_probabilities[ty]:.0%}[/]"
                name = f"[bold blue]{ty}[/]"
                self.warn(f"Couldn't assign type to [bold blue]{cluster.id}[/] (maybe {name}, {prob} confidence)")

        return clusters_new


class ClusterLoaderMixin(Command):

    clusters: str

    def _load_clusters(self) -> "ClusterTable":
        from ...model import ClusterTable

        try:
            # get filesize and unit
            input_size = os.stat(self.clusters).st_size
            total, scale, unit = ProgressReader.scale_size(input_size)
            task = self.progress.add_task("Loading clusters", total=total, unit=unit, precision=".1f")
            # load clusters
            self.info("Loading", "clusters table from file", repr(self.clusters))
            with zopen(self.clusters) as file:
                with ProgressReader(file, self.progress, task, scale) as in_:
                    return ClusterTable.load(in_)   # type: ignore
        except FileNotFoundError as err:
            self.error("Could not find clusters file:", repr(self.clusters))
            raise CommandExit(err.errno) from err

    def _label_genes(self, genes: List["Gene"], clusters: "ClusterTable") -> List["Gene"]:
        cluster_by_seq = collections.defaultdict(list)
        for i in range(len(clusters)):
            cluster_by_seq[clusters.sequence_id[i]].append((clusters.start[i], clusters.end[i]))

        gene_count = len(genes)
        unit = "gene" if gene_count == 1 else "genes"
        task = self.progress.add_task("Labelling genes", total=gene_count, unit=unit, precision="")

        self.info("Labelling", "genes belonging to clusters")
        labelled_genes = []
        for seq_id, seq_genes in itertools.groupby(genes, key=operator.attrgetter("source.id")):
            for gene in seq_genes:
                if any(
                    cluster_start <= gene.start and gene.end <= cluster_end
                    for (cluster_start, cluster_end) in cluster_by_seq[seq_id]
                ):
                    gene = gene.with_probability(1)
                else:
                    gene = gene.with_probability(0)
                labelled_genes.append(gene)
                self.progress.update(task_id=task, advance=1)

        return labelled_genes

    def _extract_clusters(self, genes: List["Gene"], clusters: "ClusterTable") -> List["Cluster"]:
        from ...model import Cluster, ClusterType

        cluster_types = {}
        cluster_by_seq = collections.defaultdict(list)
        for i in range(len(clusters)):
            seq_id = clusters.sequence_id[i]
            cluster_id = clusters.cluster_id[i]
            cluster_by_seq[seq_id].append((
                clusters.start[i],
                clusters.end[i],
                clusters.cluster_id[i],
            ))
            if not "type" in clusters.data.columns:
                cluster_types[cluster_id] = None
            elif clusters.type[i] == "Unknown" or clusters.type[i] is None:
                cluster_types[cluster_id] = ClusterType()
            else:
                cluster_types[cluster_id] = ClusterType(*clusters.type[i].split(";"))

        self.info("Extracting", "genes belonging to clusters")
        genes_by_cluster = collections.defaultdict(list)
        for seq_id, seq_genes in itertools.groupby(genes, key=operator.attrgetter("source.id")):
            for gene in seq_genes:
                for cluster_start, cluster_end, cluster_id in cluster_by_seq[seq_id]:
                    if cluster_start <= gene.start and gene.end <= cluster_end:
                        genes_by_cluster[cluster_id].append(gene)

        return [
            Cluster(cluster_id, genes_by_cluster[cluster_id], cluster_types[cluster_id])
            for cluster_id in sorted(clusters.cluster_id)
            if genes_by_cluster[cluster_id]
        ]


class CompositionWriterMixin(Command):

    output_dir: str

    def _save_domain_compositions(self, all_possible: List[str], clusters: List["Cluster"]) -> None:
        import numpy
        import scipy.sparse

        self.info("Saving", "training matrix labels for type classifier")
        with open(os.path.join(self.output_dir, "domains.tsv"), "w") as out:
            out.writelines(f"{domain}\n" for domain in all_possible)
        with open(os.path.join(self.output_dir, "types.tsv"), "w") as out:
            writer = csv.writer(out, dialect="excel-tab")
            for cluster in clusters:
                types = ";".join(sorted(cluster.type.names))
                writer.writerow([cluster.id, types])

        self.info("Building", "new domain composition matrix")
        comp = numpy.array([
            c.domain_composition(all_possible)
            for c in clusters
        ])

        comp_out = os.path.join(self.output_dir, "compositions.npz")
        self.info("Saving", "new domain composition matrix to file", repr(comp_out))
        scipy.sparse.save_npz(comp_out, scipy.sparse.coo_matrix(comp))

