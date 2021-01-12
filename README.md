![](static/gecco.png)

# Hi, I'm GECCO!


## Requirements

* [Python](https://www.python.org/downloads/) 3.6 or higher


## Installing GECCO

### Development version

To install GECCO from the EMBL git server, run:
```console
$ pip install git+https://git.embl.de/grp-zeller/GECCO/
```

Note that this command can take a long time to complete as it need to download
around 250MB of data from the EBI FTP server. Once the install is finished, a
`gecco` command should be available in your path.


## Running GECCO

Once `gecco.py` is available in your `PATH`, you can run it from everywhere by
giving it a FASTA or GenBank file with the genome you want to analyze, as well
as an output directory.

```console
$ gecco run --genome some_genome.fna -o some_output_dir
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
