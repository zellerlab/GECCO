## Hi, I'm GECCO!

### Requirements

* [Python](https://www.python.org/downloads/) 3.6 or higher
* [HMMER](http://hmmer.org/) v3.2 or higher
* [Prodigal](https://github.com/hyattpd/Prodigal) v2.6.3 or higher

Both HMMER and Prodigal can be installed through conda. They have to be in the `PATH` variable before running GECCO.


### Installation

To install GECCO, just run

```bash
git clone https://gitlab.com/astair/gecco.git
cd GECCO
python setup.py install
```

This should install all Python dependencies and download the Pfam database. Also, it should add `gecco.py` to your `PATH` variable. If that did not work for some reason, you can do this manually by typing

```bash
export PATH=$(pwd):$PATH
```


### Running GECCO

Once `gecco.py` is available in your `PATH`, you can run it from everywhere by giving it a FASTA or GenBank file with the genome you want to analyze, as well as an output directory.

```bash
gecco.py -o some_output_dir some_genome.fna
```

For more info, you could check the wiki... if there was one. So far, you're out of luck, because I haven't gotten aroud to write one yet ;)
