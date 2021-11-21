Hi, I’m GECCO!
==============

.. image:: https://raw.githubusercontent.com/zellerlab/GECCO/v0.6.2/static/gecco-square.png
   :target: https://github.com/zellerlab/GECCO/


🦎 ️Overview
---------------

GECCO (Gene Cluster prediction with Conditional Random Fields) is a fast
and scalable method for identifying putative novel Biosynthetic Gene
Clusters (BGCs) in genomic and metagenomic data using Conditional Random
Fields (CRFs).

|GitLabCI| |License| |Coverage| |Docs| |Source| |Mirror| |Changelog|
|Issues| |Preprint| |PyPI| |Bioconda| |Versions| |Wheel|

🔧 Installing GECCO
-------------------

GECCO is implemented in `Python <https://www.python.org/>`__, and
supports `all versions <https://endoflife.date/python>`__ from Python
3.6. It requires additional libraries that can be installed directly
from `PyPI <https://pypi.org>`__, the Python Package Index.

Use `pip <https://pip.pypa.io/en/stable/>`__ to install GECCO on
your machine::

   $ pip install gecco-tool

If you’d rather use `Conda <https://conda.io>`__, a package is available
in the `bioconda <https://bioconda.github.io/>`__ channel. You can
install with::

   $ conda install -c bioconda gecco

This will install GECCO, its dependencies, and the data needed to run
predictions. This requires around 100MB of data to be downloaded, so it
could take some time depending on your Internet connection. Once done,
you will have a ``gecco`` command available in your $PATH.

*Note that GECCO uses* `HMMER3 <http://hmmer.org/>`__, *which can
only run on PowerPC and recent x86-64 machines running a POSIX operating
system. Therefore, Linux and OSX are supported platforms, but GECCO will
not be able to run on Windows.*

🧬 Running GECCO
-----------------

Once ``gecco`` is installed, you can run it from the terminal by giving
it a FASTA or GenBank file with the genomic sequence you want to
analyze, as well as an output directory::

   $ gecco run --genome some_genome.fna -o some_output_dir

Additional parameters of interest are:

-  ``--jobs``, which controls the number of threads that will be spawned
   by GECCO whenever a step can be parallelized. The default, *0*, will
   autodetect the number of CPUs on the machine using
   `os.cpu_count <https://docs.python.org/3/library/os.html#os.cpu_count>`__.
-  ``--cds``, controlling the minimum number of consecutive genes a BGC
   region must have to be detected by GECCO (default is 3).
-  ``--threshold``, controlling the minimum probability for a gene to be
   considered part of a BGC region. Using a lower number will increase
   the number (and possibly length) of predictions, but reduce accuracy.

🔖 Reference
-------------

GECCO can be cited using the following preprint:

   **Accurate de novo identification of biosynthetic gene clusters with
   GECCO**. Laura M Carroll, Martin Larralde, Jonas Simon Fleck, Ruby
   Ponnudurai, Alessio Milanese, Elisa Cappio Barazzone, Georg Zeller.
   bioRxiv 2021.05.03.442509;
   `doi:10.1101/2021.05.03.442509 <https://doi.org/10.1101/2021.05.03.442509>`__

💭 Feedback
------------

⚠️ Issue Tracker
~~~~~~~~~~~~~~~~

Found a bug ? Have an enhancement request ? Head over to the `GitHub
issue tracker <https://github.com/zellerlab/GECCO/issues>`__ if you need
to report or ask something. If you are filing in on a bug, please
include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.

🏗️ Contributing
~~~~~~~~~~~~~~~~

Contributions are more than welcome! See
`CONTRIBUTING.md <https://github.com/zellerlab/GECCO/blob/master/CONTRIBUTING.md>`__
for more details.

⚖️ License
----------

This software is provided under the `GNU General Public License v3.0 or
later <https://choosealicense.com/licenses/gpl-3.0/>`__. GECCO is
developped by the `Zeller
Team <https://www.embl.de/research/units/scb/zeller/index.html>`__ at
the `European Molecular Biology Laboratory <https://www.embl.de/>`__ in
Heidelberg.

.. |GitLabCI| image:: https://img.shields.io/gitlab/pipeline/grp-zeller/GECCO/master?gitlab_url=https%3A%2F%2Fgit.embl.de&style=flat-square&maxAge=600
   :target: https://git.embl.de/grp-zeller/GECCO/-/pipelines/
.. |License| image:: https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400
   :target: https://choosealicense.com/licenses/gpl-3.0/
.. |Coverage| image:: https://img.shields.io/codecov/c/gh/zellerlab/GECCO?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/zellerlab/GECCO/
.. |Docs| image:: https://img.shields.io/badge/docs-gecco.embl.de-green.svg?maxAge=2678400&style=flat-square
   :target: https://gecco.embl.de
.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/zellerlab/GECCO/
.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/grp-zeller/GECCO/
.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/zellerlab/GECCO/blob/master/CHANGELOG.md
.. |Issues| image:: https://img.shields.io/github/issues/zellerlab/GECCO.svg?style=flat-square&maxAge=600
   :target: https://github.com/zellerlab/GECCO/issues
.. |Preprint| image:: https://img.shields.io/badge/preprint-bioRxiv-darkblue?style=flat-square&maxAge=2678400
   :target: https://www.biorxiv.org/content/10.1101/2021.05.03.442509v1
.. |PyPI| image:: https://img.shields.io/pypi/v/gecco-tool.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/gecco-tool
.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/gecco?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/gecco
.. |Versions| image:: https://img.shields.io/pypi/pyversions/gecco-tool.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/gecco-tool/#files
.. |Wheel| image:: https://img.shields.io/pypi/wheel/gecco-tool?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/gecco-tool/#files
