# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
[Unreleased]: https://github.com/zellerlab/GECCO/compare/v0.10.1...master

## [v0.10.1] - 2025-11-26
[v0.10.1]: https://github.com/zellerlab/GECCO/compare/v0.10.0...v0.10.1

### Changed
- Remove unneeded `psutil` package dependency.

### Fixed
- Missing `--output` flag in `gecco convert` command ([#19](https://github.com/zellerlab/GECCO/issues/19)).


## [v0.10.0] - 2025-10-17
[v0.10.0]: https://github.com/zellerlab/GECCO/compare/v0.9.10...v0.10.0

### Added
- Support for changing the CLI argument defaults from a call to the `gecco.cli.main` function.

### Changed
- Allow training with sequences only containing negative samples (with a warning).
- Add an error in `gecco convert` when given input is not a folder.
- Update `numpy`, `polars`, `pyhmmer` and `psutil` to latest versions.
- Replace `docopt` with `argparse` parser for CLI argument parsing.
- Make `embedded_hmms` generic over the module name.
- Make `ClusterCRF.trained` accept any traversable.
- Allow changing the default `ClusterCRF` and `TypeClassifier` classes used in `gecco.cli`.

### Fixed
- `gecco.cli` and `gecco.crf` submodules not being included in wheel distribution.
- `gecco convert` issue when given a relative path.
- Unused `tqdm` import in `gecco.crf` module.
- Unused build dependencies from `setup.cfg` and `pyproject.toml`.
- `rich` optional import logic in `setup.py`.


## [v0.9.10] - 2024-02-27
[v0.9.10]: https://github.com/zellerlab/GECCO/compare/v0.9.9...v0.9.10

### Fixed
- Progress reading display when reading from compressed files.
- Change labeling routine to use broad overlaps when annotating genes with cluster tables ([#15](https://github.com/zellerlab/GECCO/pull/15)).

### Changed
- Bump supported `polars` dependency to `v0.20`.
- Bump supported `statsmodels` dependency to `v0.14`.
- Report identifier of sequences with uni-valued labels when training.


## [v0.9.9] - 2023-11-23
[v0.9.9]: https://github.com/zellerlab/GECCO/compare/v0.9.8...v0.9.9

### Added
- Support for `gzip`, `bzip2`, `lz4` and `xz`-compressed input files.

### Fixed
- Outdated use of `pandas` API in `gecco cv` command.

### Changed
- Bump `pyhmmer` dependency to `v0.10.0`.
- Bump `pyrodigal` dependency to `v3.0.0`.
- Make `gecco cv` output a gene table with a ground truth column.


## [v0.9.8] - 2023-06-09
[v0.9.8]: https://github.com/zellerlab/GECCO/compare/v0.9.7...v0.9.8

### Fixed
- `ClusterTable.from_clusters` extracting cluster IDs in the wrong column.
- Deprecation warnings in `polars.read_csv` and `polars.write_csv` with recent `polars` versions.
- Deprecation warnings in `importlib_resources` with recent Python versions.


## [v0.9.7] - 2023-05-26
[v0.9.7]: https://github.com/zellerlab/GECCO/compare/v0.9.6...v0.9.7

### Added
- Command line option to annotate proteins using bitscore cutoffs from HMMs.
- Command line option to disentangle overlapping domains after HMM annotation.

### Changed
- Bump `pyhmmer` dependency to `v0.8.0`.
- Bump `pyrodigal` dependency to `v2.1.0`.
- Rewrite `gecco.model` to use `polars` for managing tabular data.
- Replace `pandas` dependencies with `polars`
- Update `gecco run` to skip type classification for tasks without an assigned cluster type.

### Fixed
- `Cluster.to_seq_record` crashing when called on a cluster with `types` attribute unset.
- Progress bar resetting when performing domain annotation with multiple HMMs.

### Removed
- Support for Python 3.7.


## [v0.9.6] - 2023-01-11
[v0.9.6]: https://github.com/zellerlab/GECCO/compare/v0.9.5...v0.9.6

### Added
- Gene Ontology annotations to `gecco.interpro` local metadata.
- Reference to Gene Ontology terms and derived functions to `gecco.model.Domain` objects.
- Gene color based on predicted function in `gecco.model.Gene.to_seq_feature`.

### Fixed
- Missing `gzip` import in the CLI preventing usage of gzip-compressed inputs.
- Invalid coordinates of domains found in reverse-strand genes.
- Detection of entry points with `importlib.metadata` on older Python versions.

### Changed
- `bgc_id` columns of cluster tables are renamed `cluster_id`.
- `gecco.model.ProductType` is renamed to `gecco.model.ClusterType`.
- Bumped `pyrodigal` dependency to `v2.0`.
- Bumped `pyhmmer` dependency to `v0.7`.


## [v0.9.5] - 2022-08-10
[v0.9.5]: https://github.com/zellerlab/GECCO/compare/v0.9.4...v0.9.5

### Added
- `gecco predict` command to predict BGCs from an annotated genome.
- `Protein.with_seq` function to assign a new sequence to a protein object.

### Fixed
- Issue with antiSMASH sideload JSON file generation in `gecco run` and `gecco predict`.
- Make `gecco.orf` handle STOP codons consistently ([#9](https://github.com/zellerlab/GECCO/issues/9)).


## [v0.9.4] - 2022-05-31
[v0.9.4]: https://github.com/zellerlab/GECCO/compare/v0.9.3...v0.9.4

### Added
- `classes_` property to `TypeClassifier` to access the `classes_` attribute of the `TypeBinarizer`.
- Alternative ORF finder `CDSFinder` which simply extracts CDS features from input sequences ([#8](https://github.com/zellerlab/GECCO/issues/8)).
- Support for annotating domains with "exclusive" HMMs to annotate genes with *at most* one HMM from the library.

### Changed
- `ProductType` is not restricted to MIBiG types anymore and can support any string as a base type identifier.
- `PyrodigalFinder` now uses `multiprocessing.pool.ThreadPool` instead of custom thread code thanks to `OrfFinder.find_genes` reentrancy introduced in Pyrodigal `v1.0`. 
- `PyrodigalFinder` can now be used in single / non-meta mode from the API.
- BUmped minimum `rich` version to `12.3` to use `None` total in progress bars when the size of an HMM library is unknown.

### Fixed
- Broken MyPy type annotations in the `gecco.model` and `gecco.cli` modules.


## [v0.9.3] - 2022-05-13
[v0.9.3]: https://github.com/zellerlab/GECCO/compare/v0.9.2...v0.9.3

### Changed
- `--format` flag of `gecco annotate` and `gecco run` CLI commands is now made lowercase before giving value to `Bio.SeqIO`.

### Fixed
- Genes with duplicate IDs being silently ignored in `HMMER.run`.


## [v0.9.2] - 2022-04-11
[v0.9.2]: https://github.com/zellerlab/GECCO/compare/v0.9.1...v0.9.2

### Added
- Padding of short sequences with empty genes when predicting probabilities in `ClusterCRF`.

## [v0.9.1] - 2022-04-05
[v0.9.1]: https://github.com/zellerlab/GECCO/compare/v0.9.1-alpha4...v0.9.1

### Changed
- Make the `genes.tsv` and `features.tsv` table contain all genes even when they come from a contig too short to be processed by the CRF sliding window.
- Replaced the `--force-clusters-tsv` flag with a `--force-tsv` flag to force writing TSV tables even when no genes or clusters were found in `gecco run` or `gecco annotate`.

## [v0.9.1-alpha4] - 2022-03-31
[v0.9.1-alpha4]: https://github.com/zellerlab/GECCO/compare/v0.9.1-alpha3...v0.9.1-alpha4

Retrain internal model with:
```
$ python -m gecco -vv train --c1 0.4 --c2 0 --select 0.25 --window-size 20 \
         -f mibig-2.0.proG2.Pfam-v35.0.features.tsv \
         -c mibig-2.0.proG2.clusters.tsv \
         -g GECCO-data/data/embeddings/mibig-2.0.proG2.genes.tsv \
         -o models/v0.9.1-alpha4
```

## [v0.9.1-alpha3] - 2022-03-23
[v0.9.1-alpha3]: https://github.com/zellerlab/GECCO/compare/v0.9.1-alpha2...v0.9.1-alpha3

### Added
- `gecco.model.GeneTable` class to store gene coordinates independently of protein domains.

### Changed
- Refactored implementation of `load` and `dump` methods for `Table` classes into a dedicated base class.
- `gecco run` and `gecco annotate` now output a gene table in addition to the feature and cluster tables.
- `gecco train` expects a gene table instead of a GFF file for the gene coordinates.

## [v0.9.1-alpha2] - 2022-03-23
[v0.9.1-alpha2]: https://github.com/zellerlab/GECCO/compare/v0.9.1-alpha1...v0.9.1-alpha2

### Fixed
- `TypeClassifier.trained` not being able to read unknown types from type tables.

## [v0.9.1-alpha1] - 2022-03-20
[v0.9.1-alpha1]: https://github.com/zellerlab/GECCO/compare/v0.8.10...v0.9.1-alpha1
Candidate release with support for a sliding window in the CRF prediction algorithm.

## [v0.8.10] - 2022-02-23
[v0.8.10]: https://github.com/zellerlab/GECCO/compare/v0.8.9...v0.8.10
### Fixed
- `--antismash-sideload` flag of `gecco run` causing command to crash.

## [v0.8.9] - 2022-02-22
[v0.8.9]: https://github.com/zellerlab/GECCO/compare/v0.8.8...v0.8.9
### Removed
- Prediction and support for the *Other* biosynthetic type of MIBiG clusters.

## [v0.8.8] - 2022-02-21
[v0.8.8]: https://github.com/zellerlab/GECCO/compare/v0.8.7...v0.8.8
### Fixed
- `ClusterRefiner` filtering method for edge genes not working as intended.
- `gecco run` and `gecco annotate` commands crashing on missing input files instead of nicely rendering the error.

## [v0.8.7] - 2022-02-18
[v0.8.7]: https://github.com/zellerlab/GECCO/compare/v0.8.6...v0.8.7
### Fixed
- `interpro.json` metadata file not being included in distribution files.
- Missing docstring for `Protein.with_domains` method.
### Changed
- Bump minimum `scikit-learn` version to `v1.0` for Python3.7+.

## [v0.8.6] - 2022-02-17 - YANKED
[v0.8.6]: https://github.com/zellerlab/GECCO/compare/v0.8.5...v0.8.6
### Added
- CLI flag for enabling region masking for contigs processed by Prodigal.
- CLI flag for controlling region distance used for edge distance filtering.
### Changed
- `gecco.model.Gene` and `gecco.model.Protein` are now immutable data classes.
- Bump minimum `pyrodigal` version to `v0.6.4` to use region masking.
- Implement filtering for extracted clusters based on distance to the contig edge.
- Store InterPro metadata file uncompressed for version-control integration.
### Fixed
- Mark `BGC0000930` as `Terpene` in the type classifier data.
- Progress bar messages are now in consistent format.

## [v0.8.5] - 2021-11-21
[v0.8.5]: https://github.com/zellerlab/GECCO/compare/v0.8.4...v0.8.5
### Added
- Minimal compatibility support for running GECCO inside of Galaxy workflows.

## [v0.8.4] - 2021-09-26
[v0.8.4]: https://github.com/zellerlab/GECCO/compare/v0.8.3-post1...v0.8.4
### Fixed
- `gecco convert gbk --format bigslice` failing to run because of outdated code ([#5](https://github.com/zellerlab/GECCO/issues/5)).
- `gecco convert gbk --format bigslice` not creating files with names conforming to BiG-SLiCE expected input.
### Changed
- Bump minimum `pyrodigal` version to `v0.6.2` to use platform-accelerated code if supported.

## [v0.8.3-post1] - 2021-08-23
[v0.8.3-post1]: https://github.com/zellerlab/GECCO/compare/v0.8.3...v0.8.3-post1
### Fixed
- Wrong default value for `--threshold` being shown in `gecco run` help message.

## [v0.8.3] - 2021-08-23
[v0.8.3]: https://github.com/zellerlab/GECCO/compare/v0.8.2...v0.8.3
### Changed
- Default probability threshold for segmentation to 0.3 (from 0.4).

## [v0.8.2] - 2021-07-31
[v0.8.2]: https://github.com/zellerlab/GECCO/compare/v0.8.1...v0.8.2
### Fixed
- `gecco run` crashing on Python 3.6 because of missing `contextlib.nullcontext` class.
### Changed
- `gecco run` and `gecco annotate` will not try to count the number of profiles when given an external HMM file with the `--hmm` flag.
- `PyHMMER.run` now reports the *p-value* of each domain in addition to the *e-value* as a `/note` qualifier.

## [v0.8.1] - 2021-07-29
[v0.8.1]: https://github.com/zellerlab/GECCO/compare/v0.8.0...v0.8.1
### Changed
- `gecco run` now filters out unneeded features before annotating, making it easier to analyze the results of a run with a custom `--model`.
### Fixed
- `gecco` reporting about using Pfam `v33.1` while actually using `v34.0` because of an outdated field in `gecco/hmmer/Pfam.ini`.
### Added
- Missing documentation for the `strand` attribute of `gecco.model.Gene`.

## [v0.8.0] - 2021-07-03
[v0.8.0]: https://github.com/zellerlab/GECCO/compare/v0.7.0...v0.8.0
### Changed
- Retrain internal model using new sequence embeddings and remove broken/duplicate BGCs from MIBiG 2.0.
- Bump minimum `pyhmmer` version to `v0.4.0` to improve exception handling.
- Bump minimum `pyrodigal` version to `v0.5.0` to fix sequence decoding on some platforms.
- Use p-values instead of e-values to filter domains obtained with HMMER.
- `gecco cv` and `gecco train` now seed the RNG with a user-defined seed before shuffling rows of training data.
### Fixed
- Extraction of BGC compositions for the type predictor while training.
- `ClusterCRF.trained` failing to open an external model.
### Added
- `Domain.pvalue` attribute to access the p-value of a domain annotation.
- Mandatory `pvalue` column to `FeatureTable` objects.
- Support for loading several feature tables in `gecco train` and `gecco cv`.
- Warnings to `ClusterCRF.fit` when selecting uninformative features.
- `--correction` flag to `gecco train` and `gecco cv`, allowing to give a multiple testing correction method when computing p-values with the Fisher Exact Tests.
### Removed
- Outdated `gecco embed` command.
- Unused `--truncate` flag from the `gecco train` CLI.
- Tigrfam domains, which is not improving performance on the new training data.

## [v0.7.0] - 2021-05-31
[v0.7.0]: https://github.com/zellerlab/GECCO/compare/v0.6.3...v0.7.0
### Added
- Support for writing an AntiSMASH sideload JSON file after a `gecco run` workflow.
- Code for converting GenBank files in BiG-SLiCE compatible format with the `gecco convert` subcommand.
- Documentation about using GECCO in combination with AntiSMASH or BiG-SLiCE.
### Changed
- Minimum Biopython version to `v1.73` for compatibility with older bioinformatics tooling.
- Internal domain composition shipped in the `gecco.types` with newer composition array obtained directly from MIBiG files.
### Removed
- Outdated notice about `-vvv` verbosity level in the help message of the main `gecco` command.

## [v0.6.3] - 2021-05-10
[v0.6.3]: https://github.com/zellerlab/GECCO/compare/v0.6.2...v0.6.3
### Fixed
- HMMER annotation not properly handling inputs with multiple contigs.
- Some progress bar totals displaying as floats in the CLI.
### Changed
- `PyHMMER` now sets the `Z` and `domZ` values from the number of proteins given to the search pipeline.
- `gecco.cli` delegates imports to make CLI more responsive.
- `pkg_resources` has been replaced with `importlib.resources` and `importlib.metadata` where applicable.
- `multiprocessing.cpu_count` has been replaced with `os.cpu_count` where applicable.

## [v0.6.2] - 2021-05-04
[v0.6.2]: https://github.com/zellerlab/GECCO/compare/v0.6.1...v0.6.2
### Fixed
- `gecco cv loto` crashing because of outdated code.
### Changed
- Logging-style prompt will only display if GECCO is running with `-vv` flag.
### Added
- GECCO bioRxiv paper reference to `Cluster.to_seq_record` output record.

## [v0.6.1] - 2021-03-15
[v0.6.1]: https://github.com/zellerlab/GECCO/compare/v0.6.0...v0.6.1
### Fixed
- Progress bar not being disabled by `-q` flag in CLI.
- Fallback to using HMM name if accession is not available in `PyHMMER`.
- Group genes by source contig and process them separately in `PyHMMER` to avoid bogus E-values.
### Added
- `psutil` dependency to get the number of physical CPU cores on the host machine.
- Support for using an arbitrary mapping of positives to negatives in `gecco embed`.
### Removed
- Unused and outdated `HMMER` and `DomainRow` classes from `gecco.hmmer`.

## [v0.6.0] - 2021-02-28
[v0.6.0]: https://github.com/zellerlab/GECCO/compare/v0.5.5...v0.6.0
### Changed
- Updated internal model with a cleaned-up version of the MIBiG-2.0
  Pfam-33.1/Tigrfam-15.0 embedding.
- Updated internal InterPro catalog.
### Fixed
- Features not being grouped together in `gecco cv` and `gecco train`
  when provided with a feature table where rows were not sorted by
  protein IDs.

## [v0.5.5] - 2021-02-28
[v0.5.5]: https://github.com/zellerlab/GECCO/compare/v0.5.4...v0.5.5
### Fixed
- `gecco cv` bug causing only the last fold to be written.

## [v0.5.4] - 2021-02-28
[v0.5.4]: https://github.com/zellerlab/GECCO/compare/v0.5.3...v0.5.4
### Changed
- Replaced `verboselogs`, `coloredlogs` and `better-exceptions` with `rich`.
### Removed
- `tqdm` training dependency.
### Added
- `gecco annotate` command to produce a feature table from a genomic file.
- `gecco embed` to embed BGCs into non-BGC regions using feature tables.

## [v0.5.3] - 2021-02-21
[v0.5.3]: https://github.com/zellerlab/GECCO/compare/v0.5.2...v0.5.3
### Fixed
- Coordinates of genes in output GenBank files.
- Potential issue with the number of CPUs in `PyHMMER.run`.
### Changed
- Bump required `pyrodigal` version to `v0.4.2` to fix buffer overflow.

## [v0.5.2] - 2021-01-29
[v0.5.2]: https://github.com/zellerlab/GECCO/compare/v0.5.1...v0.5.2
### Added
- Support for downloading HMM files directly from GitHub releases assets.
- Validation of filtered HMMs with MD5 checksum.
### Fixed
- Invalid coordinates of protein domains in GenBank output files.
- `gecco.interpro` module not being added to wheel distribution.
### Changed
- Bump required `pyhmmer` version to `v0.2.1`.

## [v0.5.1] - 2021-01-15
[v0.5.1]: https://github.com/zellerlab/GECCO/compare/v0.5.0...v0.5.1
### Fixed
- `--hmm` flag being ignored in in `gecco run` command.
- `PyHMMER` using HMM names instead of accessions, causing issues with Pfam HMMs.

## [v0.5.0] - 2021-01-11
[v0.5.0]: https://github.com/zellerlab/GECCO/compare/v0.4.5...v0.5.0
### Added
- Explicit support for Python 3.9.
### Changed
- [`pyhmmer`](https://pypi.org/project/pyhmmer) is used to annotate protein sequences instead of HMMER3 binary `hmmsearch`.
- HMM files are stored in binary format to speedup parsing and reduce storage size.
- `tqdm` is now a *training*-only dependency.
- `gecco cv` now requires *training* dependencies.

## [v0.4.5] - 2020-11-23
[v0.4.5]: https://github.com/zellerlab/GECCO/compare/v0.4.4...v0.4.5
### Added
- Additional `fold` column to cross-validation table output.
### Changed
- Use sequence ID instead of protein ID to extract type from cluster in `gecco cv`.
- Install HMM data in pre-pressed format to make `hmmsearch` runs faster on short sequences.
- `gecco.orf` was rewritten to extract genes from input sequences in parallel.

## [v0.4.4] - 2020-09-30
[v0.4.4]: https://github.com/zellerlab/GECCO/compare/v0.4.3...v0.4.4
### Added
- `gecco cv loto` command to run LOTO cross-validation using BGC types
  for stratification.
- `header` keyword argument to `FeatureTable.dump` and `ClusterTable.dump`
  to write the table without the column header allowing to append to an
  existing table.
- `__getitem__` implementation for `FeatureTable` and `ClusterTable`
  that returns a single row or a sub-table from a table.
### Fixed
- `gecco cv` command now writes results iteratively instead of holding
  the tables for every fold in memory.
### Changed
- Bumped `pandas` training dependency to `v1.0`.

## [v0.4.3] - 2020-09-07
[v0.4.3]: https://github.com/zellerlab/GECCO/compare/v0.4.2...v0.4.3
### Fixed
- GenBank files being written with invalid `/cds` feature type.
### Changed
- Blocked installation of Biopython `v1.78` or newer as it removes `Bio.Alphabet`
  and breaks the current code.

## [v0.4.2] - 2020-08-07
[v0.4.2]: https://github.com/zellerlab/GECCO/compare/v0.4.1...v0.4.2
### Fixed
- `TypeClassifier.predict_types` using inverse type probabilities when
  given several clusters to process.

## [v0.4.1] - 2020-08-07
[v0.4.1]: https://github.com/zellerlab/GECCO/compare/v0.4.0...v0.4.1
### Fixed
- `gecco run` command crashing on input sequences not containing any genes.

## [v0.4.0] - 2020-08-06
[v0.4.0]: https://github.com/zellerlab/GECCO/compare/v0.3.0...v0.4.0
### Added
- `gecco.model.ProductType` enum to model the biosynthetic class of a BGC.
### Removed
- `pandas` interaction from internal data model.
- `ClusterCRF` code specific to cross-validation.
### Changed
- `pandas`, `fisher` and `statsmodels` dependencies are now optional.
- `gecco train` command expects a cluster table in addition to the feature
   table to know the types of the input BGCs.

## [v0.3.0] - 2020-08-03
[v0.3.0]: https://github.com/zellerlab/GECCO/compare/v0.2.2...v0.3.0
### Changed
- Replaced Nearest-Neighbours classifier with Random Forest to perform type
  prediction for candidate BGCs.
- `gecco.knn` module was renamed to implementation-agnostic name `gecco.types`.
### Fixed
- Extraction of domain composition taking a long time in `gecco train` command.
### Removed
- `--metric` argument to the `gecco run` CLI command.

## [v0.2.2] - 2020-07-31
[v0.2.2]: https://github.com/zellerlab/GECCO/compare/v0.2.1...v0.2.2
### Changed
- `Domain` and `Gene` can now carry qualifiers that are used when they
  are translated to a sequence feature.
### Added
- InterPro names, accessions, and HMMER e-value for each annotated domain
  in GenBank output files.

## [v0.2.1] - 2020-07-23
[v0.2.1]: https://github.com/zellerlab/GECCO/compare/v0.2.0...v0.2.1
### Fixed
- Various potential crashes in `ClusterRefiner` code.
### Removed
- Uneeded feature dictionary filtering in `ClusterCRF` for models with
  Fisher Exact Test feature selection.

## [v0.2.0] - 2020-07-23
[v0.2.0]: https://github.com/zellerlab/GECCO/compare/v0.1.1...v0.2.0
### Fixed
- `pandas` warning about unsorted columns in `gecco run`.
### Removed
- `Gene.probability` property, replaced by `Gene.maximum_probability` and
  `Gene.average_probability` properties to be explicit.
### Changed
- Internal model now uses `Pfam` and `Tigrfam` with the top 35% features
  selected with Fisher's Exact Test.
- `ClusterRefiner` now removes genes on `Cluster` edges if they do not
  contain any domain annotation.

## [v0.1.1] - 2020-07-22
[v0.1.1]: https://github.com/zellerlab/GECCO/compare/v0.1.0...v0.1.1
### Added
- `ClusterCRF.predict_probabilities` to annotate a list of `Gene`.
### Changed
- BGC probability is now stored at the `Domain` level instead of at the `Gene`
  level, independently of the feature extraction level used by the CRF.
- `ClusterKNN` will use the model path provided to `gecco run` if any.
### Docs
- Added this changelog file to document changes in the code.
- Added documentation to `gecco` submodules missing some.
- Included the `CHANGELOG.md` file to the generated docs.

## [v0.1.0] - 2020-07-17
[v0.1.0]: https://github.com/zellerlab/GECCO/compare/v0.0.1...v0.1.0
Initial release.

## [v0.0.1] - 2018-08-13
[v0.0.1]: https://github.com/zellerlab/GECCO/compare/37afb97...v0.0.1
Proof-of-concept.
