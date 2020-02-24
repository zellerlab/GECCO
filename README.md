![](static/gecco.png)

# Hi, I'm GECCO!


## Requirements

* [Python](https://www.python.org/downloads/) 3.6 or higher
* [HMMER](http://hmmer.org/) v3.2 or higher
* [Prodigal](https://github.com/hyattpd/Prodigal) v2.6.3 or higher

Both HMMER and Prodigal can be installed through [conda](https://anaconda.org/).
They have to be in the `$PATH` variable before running GECCO.


## Installation

### Development version

To install GECCO from the EMBL git server, run:
```console
$ pip install git+https://git.embl.de/grp-zeller/GECCO/
```

Note that this command can take a long time to complete as it need to download
around 250MB of data from the EBI FTP server. You will need to have writing
rights to the site folder of Python; if this is not the case, use `pip` with
the `--user` flag to install it to a local folder. Another option is to use
a virtual environment, either with `virtualenv` or with `conda`.

Once the install is finished, a `gecco` command will be available in your path
automatically.


## Running GECCO

Once `gecco.py` is available in your `PATH`, you can run it from everywhere by
giving it a FASTA or GenBank file with the genome you want to analyze, as well
as an output directory.

```console
$ gecco run --genome some_genome.fna -o some_output_dir
```


## Training GECCO

### Resources

For this, you need to get the FASTA file from MIBiG containing all the proteins
of BGCs. It can be downloaded [from there](https://mibig.secondarymetabolites.org/download).

Then you also need some bacterial genomes free of BGCs, also in FASTA format. A
common approach is to download a bunch of complete genomes from the ENA and to
remove the ones where a BGC is detected with either DeepBGC or AntiSMASH.


### Build feature tables

With your training sequences ready, first build a feature table using the
`gecco annotate` subcommand:

```console
$ gecco annotate --genome reference_genome.fna -o nobgc
$ gecco annotate --proteins mibig_prot_seqs_xx.faa -o bgc
```

Use the `--hmm` flag to give other HMMs instead of using the internal one
(PFam v31.0 only at the moment).

*This step will probably take ages, count about 5 minutes to annotate
1M amino acids with PFam.*


### Create the embedding

When both the BGC and the non BGC sequences have been annotated, merge them into
a continuous table to train the CRF on using the `gecco embed` subcommand:

```console
$ gecco embed --bgc bgc/*.features.tsv --nobgc nobgc/*.features.tsv -o merged.features.tsv
```

If the non-BGC training set is not large enough to fit all BGCs, a warning will
be thrown.


### Train the model

To train the model with default parameters, use the following command:

```console
$ gecco train -i merged.features.tsv -o model
```
