<img align="right" width="180" height="180" src="https://raw.githubusercontent.com/zellerlab/GECCO/v0.6.2/static/gecco-square.png">

# Hi, I'm GECCO!

## 🦎 ️Overview

GECCO (Gene Cluster prediction with Conditional Random Fields) is a fast and
scalable method for identifying putative novel Biosynthetic Gene Clusters (BGCs)
in genomic and metagenomic data using Conditional Random Fields (CRFs).

[![Actions](https://img.shields.io/github/actions/workflow/status/zellerlab/GECCO/test.yml?branch=master&style=flat-square&maxAge=300)](https://github.com/zellerlab/GECCO/actions/workflows/test.yml)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![Coverage](https://img.shields.io/codecov/c/gh/zellerlab/GECCO?style=flat-square&maxAge=600)]( https://codecov.io/gh/zellerlab/GECCO/)
[![Docs](https://img.shields.io/badge/docs-gecco.embl.de-green.svg?maxAge=2678400&style=flat-square)](https://gecco.embl.de)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/zellerlab/GECCO/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/grp-zeller/GECCO/)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/zellerlab/GECCO/blob/master/CHANGELOG.md)
[![Issues](https://img.shields.io/github/issues/zellerlab/GECCO.svg?style=flat-square&maxAge=600)](https://github.com/zellerlab/GECCO/issues)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-darkblue?style=flat-square&maxAge=2678400)](https://www.biorxiv.org/content/10.1101/2021.05.03.442509v1)
[![PyPI](https://img.shields.io/pypi/v/gecco-tool.svg?style=flat-square&maxAge=3600)](https://pypi.python.org/pypi/gecco-tool)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/gecco?style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/gecco)
[![Galaxy](https://img.shields.io/badge/Galaxy-GECCO-darkblue?style=flat-square&maxAge=3600)](https://toolshed.g2.bx.psu.edu/repository?repository_id=c29bc911b3fc5f8c)
[![Versions](https://img.shields.io/pypi/pyversions/gecco-tool.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/gecco-tool/#files)
[![Wheel](https://img.shields.io/pypi/wheel/gecco-tool?style=flat-square&maxAge=3600)](https://pypi.org/project/gecco-tool/#files)


## 🔧 Installing GECCO

GECCO is implemented in [Python](https://www.python.org/), and supports [all
versions](https://endoflife.date/python) from Python 3.7. It requires
additional libraries that can be installed directly from
[PyPI](https://pypi.org), the Python Package Index.

Use [`pip`](https://pip.pypa.io/en/stable/) to install GECCO on your
machine:
```console
$ pip install gecco-tool
```

If you'd rather use [Conda](https://conda.io), a package is available
in the [`bioconda`](https://bioconda.github.io/) channel. You can install
with:
```console
$ conda install -c bioconda gecco
```

This will install GECCO, its dependencies, and the data needed to run
predictions. This requires around 40MB of data to be downloaded, so
it could take some time depending on your Internet connection. Once done,
you will have a ``gecco`` command available in your $PATH.

*Note that GECCO uses [HMMER3](http://hmmer.org/), which can only run
on PowerPC and recent x86-64 machines running a POSIX operating system.
Therefore, GECCO will work on Linux and OSX, but not on Windows.*


## 🧬 Running GECCO

Once `gecco` is installed, you can run it from the terminal by giving it a
FASTA or GenBank file with the genomic sequence you want to analyze, as
well as an output directory:

```console
$ gecco run --genome some_genome.fna -o some_output_dir
```

Additional parameters of interest are:

- `--jobs`, which controls the number of threads that will be spawned by
  GECCO whenever a step can be parallelized. The default, *0*, will
  autodetect the number of CPUs on the machine using
  [`os.cpu_count`](https://docs.python.org/3/library/os.html#os.cpu_count).
- `--cds`, controlling the minimum number of consecutive genes a BGC region
  must have to be detected by GECCO. The default is *3*.
- `--threshold`, controlling the minimum probability for a gene to be
  considered part of a BGC region. Using a lower number will increase the
  number (and possibly length) of predictions, but reduce accuracy. The
  default of *0.8* was selected to optimize precision/recall on a test set
  of 364 BGCs from [MIBiG 2.0](https://mibig.secondarymetabolites.org/).
- `--cds-feature`, which can be supplied a feature name to extract genes
  if the input file already contains gene annotations instead of predicting
  genes with [Pyrodigal](https://pyrodigal.readthedocs.io). A common value
  for records downloaded from GenBank is `--cds-feature CDS`.

## 🔎 Results

GECCO will create the following files:

- `{genome}.genes.tsv`: The *genes* file, containing the genes extracted
  or predicted from the input file, and per-gene BGC probabilities
  predicted by the CRF.
- `{genome}.features.tsv`: The *features* file, containing the identified
  domains in the input sequences, in tabular format.
- `{genome}.clusters.tsv`: If any were found, a *clusters* file, containing
  the coordinates of the predicted clusters along their putative biosynthetic
  type, in tabular format.
- `{genome}_cluster_{N}.gbk`: If any were found, a GenBank file per cluster,
  containing the cluster sequence annotated with its member proteins and domains.

GECCO can also convert results to other formats that may be more convenient
depending on the downstream usage. GECCO can convert results into:

- GFF3 format so they can be loaded into a genomic viewer 
  (`gecco convert clusters --format gff`).
- GenBank files with antiSMASH-style features so they can be loaded into 
  [BiG-SLiCE](https://github.com/medema-group/bigslice) for further analysis
  (`gecco convert gbk --format bigslice`).
- FASTA files with the sequences of all the predicted BGCs (`gecco convert gbk --format fna`)
  or with the sequences of all their proteins (`gecco convert gbk --format faa`).

To get a more visual way of exploring of the predictions, you
can open the GenBank files in a genome editing software like [UGENE](http://ugene.net/).
You can otherwise load the results into an AntiSMASH report: check the
[Integrations](https://gecco.embl.de/integrations.html#antismash) page of the
documentation for a step-by-step guide. 


## 🔖 Reference

GECCO can be cited using the following preprint:

> **Accurate de novo identification of biosynthetic gene clusters with GECCO**.
> Laura M Carroll, Martin Larralde, Jonas Simon Fleck, Ruby Ponnudurai, Alessio Milanese, Elisa Cappio Barazzone, Georg Zeller.
> bioRxiv 2021.05.03.442509; [doi:10.1101/2021.05.03.442509](https://doi.org/10.1101/2021.05.03.442509)


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/zellerlab/GECCO/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See [`CONTRIBUTING.md`](https://github.com/zellerlab/GECCO/blob/master/CONTRIBUTING.md)
for more details.

## ⚖️ License

This software is provided under the [GNU General Public License v3.0 *or later*](https://choosealicense.com/licenses/gpl-3.0/). GECCO is developped by the [Zeller Team](https://www.embl.de/research/units/scb/zeller/index.html)
at the [European Molecular Biology Laboratory](https://www.embl.de/) in Heidelberg.
