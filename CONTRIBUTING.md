# Contributing

For bug fixes or new features, please file an issue before submitting a pull request.
If the change is not trivial, it may be best to wait for feedback.

## Running tests

Tests are written as usual Python unit tests with the `unittest` module of the
Python standard library. Running them can be done as follow:

```console
$ python -m unittest discover -vv
```


## Managing versions

As it is a good project management practice, we should follow
[semantic versioning](https://semver.org/), so remember the following:

* As long as the model predicts the same thing, retraining/updating the model
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
version of Pfam), simply edit the INI file for that HMM in the
``gecco/hmmer`` folder.

Then clean and rebuild data files to download the latest version of
the HMMs:

```console
$ python setup.py clean build_data --inplace
```


### Upgrading the internal CRF model

After having trained a new version of the model, run the following
command to update the internal GECCO model as well as the hash signature file:

```console
$ python setup.py update_model --model <path_to_new_crf.model>
```
