# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
[Unreleased]: https://git.embl.de/grp-zeller/GECCO/compare/v0.5.5...master

## [v0.5.5] - 2021-02-28
[v0.5.5]: https://git.embl.de/grp-zeller/GECCO/compare/v0.5.4...v0.5.5
### Fixed
- `gecco cv` bug causing only the last fold to be written.

## [v0.5.4] - 2021-02-28
[v0.5.4]: https://git.embl.de/grp-zeller/GECCO/compare/v0.5.3...v0.5.4
### Changed
- Replaced `verboselogs`, `coloredlogs` and `better-exceptions` with `rich`.
### Removed
- `tqdm` training dependency.
### Added
- `gecco annotate` command to produce a feature table from a genomic file.
- `gecco embed` to embed BGCs into non-BGC regions using feature tables.

## [v0.5.3] - 2021-02-21
[v0.5.3]: https://git.embl.de/grp-zeller/GECCO/compare/v0.5.2...v0.5.3
### Fixed
- Coordinates of genes in output GenBank files.
- Potential issue with the number of CPUs in `PyHMMER.run`.
### Changed
- Bump required `pyrodigal` version to `v0.4.2` to fix buffer overflow.

## [v0.5.2] - 2021-01-29
[v0.5.2]: https://git.embl.de/grp-zeller/GECCO/compare/v0.5.1...v0.5.2
### Added
- Support for downloading HMM files directly from GitHub releases assets.
- Validation of filtered HMMs with MD5 checksum.
### Fixed
- Invalid coordinates of protein domains in GenBank output files.
- `gecco.interpro` module not being added to wheel distribution.
### Changed
- Bump required `pyhmmer` version to `v0.2.1`.

## [v0.5.1] - 2021-01-15
[v0.5.1]: https://git.embl.de/grp-zeller/GECCO/compare/v0.5.0...v0.5.1
### Fixed
- `--hmm` flag being ignored in in `gecco run` command.
- `PyHMMER` using HMM names instead of accessions, causing issues with Pfam HMMs.

## [v0.5.0] - 2021-01-11
[v0.5.0]: https://git.embl.de/grp-zeller/GECCO/compare/v0.4.5...v0.5.0
### Added
- Explicit support for Python 3.9.
### Changed
- [`pyhmmer`](https://pypi.org/project/pyhmmer) is used to annotate protein sequences instead of HMMER3 binary `hmmsearch`.
- HMM files are stored in binary format to speedup parsing and reduce storage size.
- `tqdm` is now a *training*-only dependency.
- `gecco cv` now requires *training* dependencies.

## [v0.4.5] - 2020-11-23
[v0.4.5]: https://git.embl.de/grp-zeller/GECCO/compare/v0.4.4...v0.4.5
### Added
- Additional `fold` column to cross-validation table output.
### Changed
- Use sequence ID instead of protein ID to extract type from cluster in `gecco cv`.
- Install HMM data in pre-pressed format to make `hmmsearch` runs faster on short sequences.
- `gecco.orf` was rewritten to extract genes from input sequences in parallel.

## [v0.4.4] - 2020-09-30
[v0.4.4]: https://git.embl.de/grp-zeller/GECCO/compare/v0.4.3...v0.4.4
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
[v0.4.3]: https://git.embl.de/grp-zeller/GECCO/compare/v0.4.2...v0.4.3
### Fixed
- GenBank files being written with invalid `/cds` feature type.
### Changed
- Blocked installation of Biopython `v1.78` or newer as it removes `Bio.Alphabet`
  and breaks the current code.

## [v0.4.2] - 2020-08-07
[v0.4.2]: https://git.embl.de/grp-zeller/GECCO/compare/v0.4.1...v0.4.2
### Fixed
- `TypeClassifier.predict_types` using inverse type probabilities when
  given several clusters to process.

## [v0.4.1] - 2020-08-07
[v0.4.1]: https://git.embl.de/grp-zeller/GECCO/compare/v0.4.0...v0.4.1
### Fixed
- `gecco run` command crashing on input sequences not containing any genes.

## [v0.4.0] - 2020-08-06
[v0.4.0]: https://git.embl.de/grp-zeller/GECCO/compare/v0.3.0...v0.4.0
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
[v0.3.0]: https://git.embl.de/grp-zeller/GECCO/compare/v0.2.2...v0.3.0
### Changed
- Replaced Nearest-Neighbours classifier with Random Forest to perform type
  prediction for candidate BGCs.
- `gecco.knn` module was renamed to implementation-agnostic name `gecco.types`.
### Fixed
- Extraction of domain composition taking a long time in `gecco train` command.
### Removed
- `--metric` argument to the `gecco run` CLI command.

## [v0.2.2] - 2020-07-31
[v0.2.2]: https://git.embl.de/grp-zeller/GECCO/compare/v0.2.1...v0.2.2
### Changed
- `Domain` and `Gene` can now carry qualifiers that are used when they
  are translated to a sequence feature.
### Added
- InterPro names, accessions, and HMMER e-value for each annotated domain
  in GenBank output files.

## [v0.2.1] - 2020-07-23
[v0.2.1]: https://git.embl.de/grp-zeller/GECCO/compare/v0.2.0...v0.2.1
### Fixed
- Various potential crashes in `ClusterRefiner` code.
### Removed
- Uneeded feature dictionary filtering in `ClusterCRF` for models with
  Fisher Exact Test feature selection.

## [v0.2.0] - 2020-07-23
[v0.2.0]: https://git.embl.de/grp-zeller/GECCO/compare/v0.1.1...v0.2.0
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
[v0.1.1]: https://git.embl.de/grp-zeller/GECCO/compare/v0.1.0...v0.1.1
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
[v0.1.0]: https://git.embl.de/grp-zeller/GECCO/compare/v0.0.1...v0.1.0
Initial release.

## [v0.0.1] - 2018-08-13
[v0.0.1]: https://git.embl.de/grp-zeller/GECCO/compare/37afb97...v0.0.1
Proof-of-concept.
