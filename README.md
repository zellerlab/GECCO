![](static/gecco.png)

# Hi, I'm GECCO!


## Requirements

* [Python](https://www.python.org/downloads/) 3.6 or higher
* [HMMER](http://hmmer.org/) v3.2 or higher
* [Prodigal](https://github.com/hyattpd/Prodigal) v2.6.3 or higher

Both HMMER and Prodigal can be installed through [conda](https://anaconda.org/).
They have to be in the `$PATH` variable before running GECCO.


## Installation

To install GECCO, just run:
```console
$ pip install https://git.embl.de/grp-zeller/GECCO/-/archive/refactor/GECCO-refactor.zip
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

```bash
gecco run --genome some_genome.fna -o some_output_dir
```

For more info, you could check the wiki... if there was one.
So far, you're out of luck because I haven't gotten aroud to write one yet ;)
