![](static/gecco.png)

# Hi, I'm GECCO!


## Requirements

* [Python](https://www.python.org/downloads/) 3.6 or higher
* [HMMER](http://hmmer.org/) v3.2 or higher

HMMER can be installed through [conda](https://anaconda.org/). It has to
be in the `$PATH` variable before running GECCO. GECCO also requires
additional Python libraries, but they are normally installed automatically
by `pip` or `conda` when installing GECCO.


## Installing GECCO

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
$ gecco annotate --mibig mibig_prot_seqs_xx.faa -o bgc
```

Use the `--hmm` flag to give other HMMs instead of using the internal one
(PFam v31.0 only at the moment). **Make sure to use the `--mibig` input flag
and not `--proteins` when annotating MIBiG sequences to ensure additional
metadata are properly extracted from the sequence id of each protein.**

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



## Housekeeping GECCO


## Versioning

As it is a good project management practice, we should follow
[semantic versioning](https://semver.org/), so remember the following:

* As long as the model predict the same thing, retraining/updating the model
  should be considered a non-breaking change, so you should bump the MINOR
  version of the program.
* Upgrading the internal HMMs could potentially change the output but won't
  break the program, they should be treated as non-breaking change, so you
  should bump the MINOR version of the program.
* If the model changes prediction (e.g. predicted classes change), then you
  should bump the MAJOR version of the program as it it a breaking change.
* Changes in the code should be treated following semver the usual way.
* Changed in the CLI should be treated as changed in the API (e.g. a new
  CLI option or a new command bumps the MINOR version, removal of an option
  bumps the MAJOR version).


### Upgrading the internal HMMs

To bump the version of the internal HMMs (for instance, to switch to a newer
version of Pfam), you will need to do the following:

- edit the `setup.py` file with the new URL to the HMM file.
- update the signature file in `gecco/data/hmms` with the MD5 checksum of the
  new file (this can be found online or computed locally with the `md5sums`
  command after having downloaded the file from a safe source).


### Upgrading the internal CRF model

After having trained a new version of the model, simply run the following command
to update the internal GECCO model as well as the hash signature file:

```console
$ python setup.py update_model --model <path_to_new_crf.model>
```
