# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
[Unreleased]: https://git.embl.de/grp-zeller/GECCO/compare/v0.2.1...master

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
