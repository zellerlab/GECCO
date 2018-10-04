## Hi, I'm GECCO!

### Requirements

* Python 3.6 or higher
* HMMER v3.2 or higher (http://hmmer.org/)
* Prodigal v2.6.3 or higher (https://github.com/hyattpd/Prodigal)

Both HMMER and Prodigal can be installed through conda. They have to be in the `PATH` variable before running GECCO.


### Installation

To install GECCO, just run

```bash
git clone https://gitlab.com/astair/gecco.git
cd GECCO
python setup.py install
```

This should install all Python dependencies and download the Pfam database. Also, this should add `gecco` to your `PATH` variable. If that did not work for some reason, you can do this manually by typing

```bash
export PATH=`pwd`:$PATH
```
