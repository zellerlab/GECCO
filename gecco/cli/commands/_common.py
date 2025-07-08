import collections
import io
import itertools
import os
import operator
import pathlib
import typing
import random
from typing import List, Optional, Iterator, Iterable, Collection, Type, Callable

import rich.progress
import numpy

from .._log import ConsoleLogger
from .._utils import guess_sequences_format
from ..._meta import zopen

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord
    from ...model import Gene
    from ...hmmer import HMM
    from ...crf import ClusterCRF


# --- Output files -------------------------------------------------------------


def make_output_directory(
    logger: ConsoleLogger, output_dir: pathlib.Path, outputs: List[str]
) -> None:
    # Make output directory
    logger.info("Using", "output folder", repr(str(output_dir)), level=1)
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as err:
        logger.error("Could not create output directory: {}", err)
        raise CommandExit(err.errno) from err

    # Check if output files already exist
    for output in outputs:
        if os.path.isfile(os.path.join(output_dir, output)):
            logger.warn("Output folder contains files that will be overwritten")
            break


def write_genes_table(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    genome: pathlib.Path,
    output_dir: pathlib.Path,
) -> None:
    from ...model import GeneTable

    base, _ = os.path.splitext(os.path.basename(genome))
    pred_out = os.path.join(output_dir, f"{base}.genes.tsv")
    logger.info("Writing", "gene table to", repr(pred_out), level=1)
    with open(pred_out, "wb") as f:
        GeneTable.from_genes(genes).dump(f)


def write_feature_table(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    genome: pathlib.Path,
    output_dir: pathlib.Path,
) -> None:
    from ...model import FeatureTable

    base, _ = os.path.splitext(os.path.basename(genome))
    pred_out = os.path.join(output_dir, f"{base}.features.tsv")
    logger.info("Writing", "feature table to", repr(pred_out), level=1)
    with open(pred_out, "wb") as f:
        FeatureTable.from_genes(genes).dump(f)


def write_cluster_table(
    logger: ConsoleLogger,
    clusters: List["Cluster"],
    *,
    genome: pathlib.Path,
    output_dir: pathlib.Path,
):
    from ...model import ClusterTable

    base, _ = os.path.splitext(os.path.basename(genome))
    cluster_out = os.path.join(output_dir, f"{base}.clusters.tsv")
    logger.info("Writing", "cluster table to", repr(cluster_out), level=1)
    with open(cluster_out, "wb") as out:
        ClusterTable.from_clusters(clusters).dump(out)


def write_clusters(
    logger: ConsoleLogger,
    clusters: List["Cluster"],
    *,
    genome: pathlib.Path,
    output_dir: pathlib.Path,
    merge: bool = False,
) -> None:
    from Bio import SeqIO

    if merge:
        base, _ = os.path.splitext(os.path.basename(genome))
        gbk_out = os.path.join(output_dir, f"{base}.clusters.gbk")
        records = (cluster.to_seq_record() for cluster in clusters)
        logger.info("Writing", f"all clusters to", repr(str(gbk_out)), level=1)
        SeqIO.write(records, gbk_out, "genbank")
    else:
        for cluster in clusters:
            gbk_out = os.path.join(output_dir, f"{cluster.id}.gbk")
            logger.info(
                "Writing",
                f"cluster [bold blue]{cluster.id}[/] to",
                repr(str(gbk_out)),
                level=1,
            )
            SeqIO.write(cluster.to_seq_record(), gbk_out, "genbank")


def write_sideload_json(
    logger: ConsoleLogger,
    clusters: List["Cluster"],
) -> None:
    raise NotImplementedError


# --- Load input (sequences) ---------------------------------------------------


def load_sequences(
    logger: ConsoleLogger,
    genome: pathlib.Path,
    *,
    format: Optional[str],
) -> Iterator["SeqRecord"]:
    from Bio import SeqIO

    # try:

    # guess format or use the one given in CLI
    if format is not None:
        format: Optional[str] = format.lower()
        logger.info("Using", "user-provided sequence format", repr(format), level=2)
    else:
        logger.info("Detecting", "sequence format from file contents", level=2)
        format = guess_sequences_format(genome)
        if format is None:
            raise RuntimeError(f"Failed to detect format of {genome!r}")
        logger.success("Detected", "format of input as", repr(format), level=2)

    # load file
    with rich.progress.Progress(console=logger.console) as progress:
        # load sequences
        n = 0
        logger.info("Loading", "sequences from genomic file", repr(genome), level=1)
        with progress.open(genome, "rb", description="Loading") as file:
            # with ProgressReader(file, self.progress, task, scale) as reader:
            with zopen(file) as decompressed:
                for record in SeqIO.parse(io.TextIOWrapper(decompressed), format):  # type: ignore
                    yield record
                    n += 1

    # except FileNotFoundError as err:
    #     self.error("Could not find input file:", repr(self.genome))
    #     raise CommandExit(err.errno) from err
    # except ValueError as err:
    #     self.error("Failed to load sequences:", err)
    #     raise CommandExit(getattr(err, "errno", 1)) from err
    # else:
    #     self.success("Found", n, "sequences", level=1)


# --- Load input (tables) ------------------------------------------------------


def load_genes(
    logger: ConsoleLogger,
    table_path: pathlib.Path,
) -> Iterator["Gene"]:
    from ...model import GeneTable

    with rich.progress.Progress(console=logger.console) as progress:
        # try:
        # load gene table
        logger.info("Loading", "genes table from file", repr(str(table_path)))
        with progress.open(table_path, "rb") as file:
            # with zopen(file) as decompressed:
            table = GeneTable.load(file)
        # count genes and yield gene objects
        n_genes = len(set(table.protein_id))
        unit = "gene" if n_genes == 1 else "genes"
        task = progress.add_task("Building genes", total=n_genes)
        yield from progress.track(table.to_genes(), task_id=task)
        # except OSError as err:
        # self.error("Fail to parse genes coordinates: {}", err)
        # raise CommandExit(err.errno) from err


def load_features(
    logger,
    table_paths: Iterable[pathlib.Path],
) -> "FeatureTable":
    from ...model import FeatureTable

    features = FeatureTable()
    with rich.progress.Progress(console=logger.console) as progress:
        for filename in progress.track(table_paths):
            # try:
            # load features
            logger.info("Loading", "features table from file", repr(str(filename)))
            with progress.open(filename, "rb") as file:
                # with zopen(file) as decompressed:
                features += FeatureTable.load(file)
        # except FileNotFoundError as err:
        #     self.error("Could not find feature file:", repr(filename))
        #     raise CommandExit(err.errno) from err

    logger.success("Loaded", "a total of", len(features), "features", level=1)
    return features


def annotate_genes(
    logger: ConsoleLogger, genes: List["Gene"], features: "FeatureTable"
) -> List["Gene"]:
    from ...model import Domain

    # index genes by protein_id
    gene_index = {gene.protein.id: gene for gene in genes}
    if len(gene_index) < len(genes):
        raise ValueError("Duplicate gene names in input genes")

    # add domains from the feature table
    with rich.progress.Progress(console=logger.console) as progress:
        unit = "row" if len(features) == 1 else "rows"
        task = progress.add_task("Annotating genes", total=len(features))
        for i in progress.track(range(len(features)), task_id=task):
            # get gene by ID and check consistency
            length = features.end[i] - features.start[i]
            gene = gene_index[features.protein_id[i]]
            if gene.source.id != features.sequence_id[i]:
                raise ValueError(
                    f"Mismatched source sequence for {features.protein_id[i]!r}: {gene.source.id!r} != {features.sequence_id[i]!r}"
                )
            elif gene.end - gene.start != length:
                raise ValueError(
                    f"Mismatched gene length for {row.protein_id!r}: {gene.end - gene.start!r} != {length!r}"
                )
            elif gene.start != features.start[i]:
                raise ValueError(
                    f"Mismatched gene start for {row.protein_id!r}: {gene.start!r} != {features.start[i]!r}"
                )
            elif gene.end != features.end[i]:
                raise ValueError(
                    f"Mismatched gene end for {row.protein_id!r}: {gene.end!r} != {features.end[i]!r}"
                )
            elif gene.strand.sign != features.strand[i]:
                raise ValueError(
                    f"Mismatched gene strand for {row.protein_id!r}: {gene.strand.sign!r} != {features.strand[i]!r}"
                )
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


def assign_sources(
    logger: ConsoleLogger,
    sequences: Iterable["SeqRecord"],
    genes: List["Gene"],
    *,
    genome: pathlib.Path,
) -> Iterable["Gene"]:
    from ...model import Protein, Strand

    logger.info("Extracting", "known sequences from gene list", level=2)
    known_sequences = {gene.source.id for gene in genes}

    logger.info("Building", "index of source sequences", level=2)
    sequence_index = {
        record.id: record for record in sequences if record.id in known_sequences
    }

    # try:
    logger.info("Assigning", "source sequences to gene objects", level=2)
    for gene in genes:
        gene = gene.with_source(sequence_index[gene.source.id])
        gene_seq = gene.source.seq[gene.start - 1 : gene.end]
        if gene.strand == Strand.Reverse:
            gene_seq = gene_seq.reverse_complement()
        gene = gene.with_protein(gene.protein.with_seq(gene_seq.translate()))
        yield gene
    # except KeyError as err:
    #     logger.error("Sequence {!r} not found in {!r}", gene.source.id, str(genome))
    #     raise CommandExit(1) from err


def load_clusters(
    logger: ConsoleLogger,
    clusters: pathlib.Path,
) -> "ClusterTable":
    from ...model import ClusterTable

    # try:
    # load clusters
    logger.info("Loading", "clusters table from file", repr(str(clusters)))

    with rich.progress.Progress(console=logger.console) as progress:
        with progress.open(clusters, "rb") as file:
            # with zopen(reader) as decompressed:
            return ClusterTable.load(file)  # type: ignore
    # except FileNotFoundError as err:
    #     logger.error("Could not find clusters file:", repr(self.clusters))
    #     raise CommandExit(err.errno) from err


def label_genes(
    logger: ConsoleLogger, genes: List["Gene"], clusters: "ClusterTable"
) -> List["Gene"]:

    cluster_by_seq = collections.defaultdict(list)
    for i in range(len(clusters)):
        cluster_by_seq[clusters.sequence_id[i]].append(
            (clusters.start[i], clusters.end[i])
        )

    with rich.progress.Progress(console=logger.console) as progress:
        gene_count = len(genes)
        unit = "gene" if gene_count == 1 else "genes"
        task = progress.add_task(
            "Labelling genes", total=gene_count, unit=unit, precision=""
        )

        logger.info("Labelling", "genes belonging to clusters")
        labelled_genes = []
        for seq_id, seq_genes in itertools.groupby(
            genes, key=operator.attrgetter("source.id")
        ):
            for gene in seq_genes:
                if any(
                    cluster_start <= gene.end and gene.start <= cluster_end
                    for (cluster_start, cluster_end) in cluster_by_seq[seq_id]
                ):
                    gene = gene.with_probability(1)
                else:
                    gene = gene.with_probability(0)
                labelled_genes.append(gene)
                progress.update(task_id=task, advance=1)

    return labelled_genes


# --- Extract genes ------------------------------------------------------------


def extract_genes(
    logger,
    sequences: List["SeqRecord"],
    *,
    cds_feature: Optional[str],
    locus_tag: Optional[str],
    mask: bool,
    jobs: int,
) -> Iterable["Gene"]:
    from ...orf import PyrodigalFinder, CDSFinder

    logger.info("Extracting", "genes from input sequences", level=1)
    if cds_feature is None:
        logger.info("Using", "Pyrodigal in metagenomic mode", level=2)
        orf_finder: ORFFinder = PyrodigalFinder(metagenome=True, mask=mask, cpus=jobs)
    else:
        logger.info("Using", f"record features named {cds_feature!r}", level=2)
        orf_finder = CDSFinder(feature=cds_feature, locus_tag=locus_tag)

    unit = "contigs" if len(sequences) > 1 else "contig"

    with rich.progress.Progress(console=logger.console) as progress:
        task = progress.add_task(
            description="Finding ORFs", total=len(sequences), unit=unit, precision=""
        )

        def callback(record: "SeqRecord", found: int) -> None:
            logger.success("Found", found, "genes in record", repr(record.id), level=2)
            progress.update(task, advance=1)

        return orf_finder.find_genes(sequences, progress=callback)


# --- Annotate genes -----------------------------------------------------------


def default_hmms() -> Iterator["HMM"]:
    from ...hmmer import embedded_hmms

    return embedded_hmms()


def custom_hmms(hmm_paths: List[pathlib.Path]) -> Iterable["HMM"]:
    from ...hmmer import HMM

    for path in hmm_paths:
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
            relabel_with=r"s/([^\.]*)(\..*)?/\1/",
        )


def filter_domains(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    e_filter: Optional[float] = None,
    p_filter: Optional[float] = None,
) -> List["Gene"]:
    # Filter i-evalue and p-value if required
    if e_filter is not None:
        logger.info("Excluding", "domains with e-value over", e_filter, level=1)
        key = lambda d: d.i_evalue < e_filter
        genes = [
            gene.with_protein(
                gene.protein.with_domains(filter(key, gene.protein.domains))
            )
            for gene in genes
        ]
    if p_filter is not None:
        logger.info("Excluding", "domains with p-value over", p_filter, level=1)
        key = lambda d: d.pvalue < p_filter
        genes = [
            gene.with_protein(
                gene.protein.with_domains(filter(key, gene.protein.domains))
            )
            for gene in genes
        ]
    if p_filter is not None or e_filter is not None:
        count = sum(len(gene.protein.domains) for gene in genes)
        logger.info("Using", "remaining", count, "domains", level=1)
    return genes


def disentangle(gene: "Gene") -> "Gene":
    if len(gene.protein.domains) <= 1:
        return gene

    domains_to_keep = []
    gene_domains = gene.protein.domains.copy()
    while gene_domains:
        # get a domain and collect other overlaping domains
        domain = gene_domains.pop()
        overlaps = [
            other
            for other in gene_domains
            if other.start <= domain.end and domain.start <= other.end
        ]
        # keep only overlapping domain with the best p-value
        if not overlaps or domain.pvalue < min(d.pvalue for d in overlaps):
            domains_to_keep.append(domain)
            for other in overlaps:
                gene_domains.remove(other)

    return gene.with_protein(gene.protein.with_domains(domains_to_keep))


def disentangle_domains(logger, genes: List["Gene"]) -> List["Gene"]:
    logger.info("Disentangling", "overlapping domains in each gene", level=1)
    return list(map(_disentangle, genes))


def annotate_domains(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    hmm_paths: List[pathlib.Path],
    default_hmms: Iterable["HMM"],
    whitelist: Optional[Collection[str]] = None,
    disentangle: bool = False,
    jobs: int = 0,
    bit_cutoffs: Optional[str] = None,
    e_filter: Optional[float] = None,
    p_filter: Optional[float] = None,
) -> List["Gene"]:
    from ...hmmer import PyHMMER

    logger.info("Running", "HMMER domain annotation", level=1)

    # Run all HMMs over ORFs to annotate with protein domains
    hmms = list(custom_hmms(hmm_paths) if hmm_paths else default_hmms)
    total = None if whitelist is None else len(whitelist)

    with rich.progress.Progress(console=logger.console) as progress:

        task_hmms = progress.add_task(
            description=f"Annotating domains",
            unit="HMMs",
            total=len(hmms),
            precision="",
        )
        task_domains = progress.add_task(
            description="", total=total, unit="domains", precision=""
        )
        for hmm in progress.track(hmms, task_id=task_hmms, total=len(hmms)):
            progress.update(task_domains, description=f" {hmm.id} v{hmm.version}")
            callback = lambda h, t: progress.update(task_domains, advance=1)
            logger.info(
                "Starting",
                f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]",
                level=2,
            )
            genes = PyHMMER(hmm, jobs, whitelist).run(
                genes, progress=callback, bit_cutoffs=bit_cutoffs
            )
            logger.success(
                "Finished",
                f"annotation with [bold blue]{hmm.id} v{hmm.version}[/]",
                level=2,
            )
        progress.update(task_id=task_domains, visible=False)

    # Count number of annotated domains
    count = sum(len(gene.protein.domains) for gene in genes)
    logger.success("Found", count, "domains across all proteins", level=1)

    # Disentangle if required
    if disentangle:
        genes = _disentangle_domains(logger, genes)

    # Filter i-evalue and p-value if required
    genes = filter_domains(
        logger,
        genes,
        e_filter=e_filter,
        p_filter=p_filter,
    )

    # Sort genes
    logger.info("Sorting", "genes by coordinates", level=2)
    genes.sort(key=operator.attrgetter("source.id", "start", "end"))
    for gene in genes:
        gene.protein.domains.sort(key=operator.attrgetter("start", "end"))

    return genes


# --- Predict probabilities ----------------------------------------------------


def load_model_domains(
    logger: ConsoleLogger,
    *,
    model: Optional[pathlib.Path],
    crf_type: Type["ClusterCRF"],
) -> typing.Set[str]:
    if model is None:
        logger.info("Loading", "embedded CRF pre-trained model", level=1)
    else:
        logger.info("Loading", "CRF pre-trained model from", repr(str(model)), level=1)
    model = crf_type.trained(model)

    # extract domain list from model features
    if model.significant_features:
        domains = model.significant_features
    else:
        domains = { domain for (domain, label) in crf.model.state_features_.keys() }

    logger.success("Found", len(domains), "selected features", level=2)
    return domains


def predict_probabilities(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    model: Optional[pathlib.Path],
    pad: bool,
    crf_type: Type["ClusterCRF"],
) -> List["Gene"]:
    if model is None:
        logger.info("Loading", "embedded CRF pre-trained model", level=1)
    else:
        logger.info("Loading", "CRF pre-trained model from", repr(str(model)), level=1)
    model = crf_type.trained(model)

    logger.info("Predicting", "cluster probabilitites with the model", level=1)

    with rich.progress.Progress(console=logger.console) as progress:

        unit = "batches" if len(genes) > 1 else "batch"
        task = progress.add_task("Predicting marginals")

        def progress_callback(i: int, total: int) -> None:
            progress.update(task_id=task, completed=i, total=total)

        return model.predict_probabilities(
            genes,
            pad=pad,
            progress=progress_callback,
        )


def extract_clusters(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    threshold: float,
    postproc: str,
    cds: int,
    edge_distance: int,
) -> List["Cluster"]:
    from ...refine import ClusterRefiner

    logger.info("Extracting", "predicted clusters", level=1)
    refiner = ClusterRefiner(
        threshold=threshold, criterion=postproc, n_cds=cds, edge_distance=edge_distance
    )

    with rich.progress.Progress(console=logger.console) as progress:
        total = len({gene.source.id for gene in genes})
        # unit = "contigs" if total > 1 else "contig"
        task = progress.add_task("Extracting clusters", total=total)
        clusters: List["Cluster"] = []
        gene_groups = itertools.groupby(genes, lambda g: g.source.id)
        for _, gene_group in progress.track(gene_groups, task_id=task, total=total):
            clusters.extend(refiner.iter_clusters(list(gene_group)))

    return clusters


# --- Predict types ------------------------------------------------------------


def load_type_classifier(
    logger: ConsoleLogger,
    *,
    model: Optional[pathlib.Path],
    classifier_type: Type["TypeClassifier"],
) -> "TypeClassifier":
    if model is None:
        logger.info("Loading", "type classifier from embedded model", level=2)
    else:
        logger.info("Loading", "type classifier from", repr(str(model)), level=2)
    return classifier_type.trained(model)


def predict_types(
    logger: ConsoleLogger, clusters: List["Cluster"], *, classifier: "TypeClassifier"
) -> List["Cluster"]:
    logger.info("Predicting", "gene cluster types", level=1)
    with rich.progress.Progress(console=logger.console) as progress:
        unit = "cluster" if len(clusters) == 1 else "clusters"
        task = progress.add_task("Predicting types", total=len(clusters))
        clusters_new = []
        for cluster in progress.track(clusters, task_id=task):
            clusters_new.extend(classifier.predict_types([cluster]))
            if cluster.type:
                name = "/".join(f"[bold blue]{name}[/]" for name in cluster.type.names)
                prob = "/".join(
                    f"[bold purple]{cluster.type_probabilities[name]:.0%}[/]"
                    for name in cluster.type.names
                )
                logger.success(
                    f"Predicted type of [bold blue]{cluster.id}[/] as {name} ({prob} confidence)"
                )
            else:
                ty = max(cluster.type_probabilities, key=cluster.type_probabilities.get)  # type: ignore
                prob = f"[bold purple]{cluster.type_probabilities[ty]:.0%}[/]"
                name = f"[bold blue]{ty}[/]"
                logger.warn(
                    f"Couldn't assign type to [bold blue]{cluster.id}[/] (maybe {name}, {prob} confidence)"
                )
    return clusters_new


# --- Train model --------------------------------------------------------------


def seed_rng(logger: ConsoleLogger, seed: int) -> None:
    logger.info("Seeding", "the random number generator with seed", seed, level=2)
    random.seed(seed)
    numpy.random.seed(seed)


def fit_model(
    logger: ConsoleLogger,
    genes: List["Gene"],
    *,
    feature_type: str,
    c1: float,
    c2: float,
    window_size: int,
    window_step: int,
    shuffle: bool,
    select: Optional[float],
    correction: Optional[str],
    jobs: int = 0,
    crf_type: Type["ClusterCRF"],
) -> "ClusterCRF":
    from ...crf import ClusterCRF

    logger.info("Creating", f"the CRF in [bold blue]{feature_type}[/] mode", level=1)
    logger.info("Using", f"provided hyperparameters (C1={c1}, C2={c2})", level=1)
    if select is not None:
        logger.info(
            "Selecting",
            f"features with Fisher Exact Test (threshold={select})",
            level=1,
        )
    logger.info(
        "Iterating",
        f"over features with a sliding window (W={window_size}, step={window_step})",
        level=1,
    )
    crf = crf_type(
        feature_type,
        algorithm="lbfgs",
        window_size=window_size,
        window_step=window_step,
        c1=c1,
        c2=c2,
    )
    logger.info("Fitting", "the CRF model to the training data")
    crf.fit(
        genes, select=select, shuffle=shuffle, correction_method=correction, cpus=jobs
    )
    return crf
