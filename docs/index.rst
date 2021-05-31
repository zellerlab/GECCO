GECCO
=====

*Biosynthetic Gene Cluster prediction with Conditional Random Fields.*

|GitLabCI| |License| |Coverage| |Source| |Mirror| |Issues| |Preprint| |PyPI| |Bioconda| |Versions| |Wheel|

.. |GitLabCI| image:: https://img.shields.io/gitlab/pipeline/grp-zeller/GECCO/master?gitlab_url=https%3A%2F%2Fgit.embl.de&logo=gitlab&style=flat-square&maxAge=600
   :target: https://git.embl.de/grp-zeller/GECCO/-/pipelines

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/zellerlab/GECCO?logo=codecov&style=flat-square&maxAge=600
   :target: https://codecov.io/gh/zellerlab/GECCO/

.. |License| image:: https://img.shields.io/badge/license-GPLv3-blue.svg?logo=gnu&style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/gpl-3.0/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?logo=git&maxAge=2678400&style=flat-square
   :target: https://github.com/zellerlab/GECCO/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?logo=git&style=flat-square&maxAge=2678400
   :target: https://git.embl.de/grp-zeller/GECCO/

.. |Issues| image:: https://img.shields.io/github/issues/zellerlab/GECCO.svg?logo=github&style=flat-square&maxAge=600
   :target: https://github.com/zellerlab/GECCO/issues

.. |Preprint| image:: https://img.shields.io/badge/preprint-bioRxiv-darkblue?style=flat-square&maxAge=2678400&logo=arxiv
   :target: https://www.biorxiv.org/content/10.1101/2021.05.03.442509v1

.. |PyPI| image:: https://img.shields.io/pypi/v/gecco-tool.svg?style=flat-square&maxAge=3600&logo=pypi
   :target: https://pypi.python.org/pypi/gecco-tool

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/gecco?logo=anaconda&style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/gecco

.. |Versions| image:: https://img.shields.io/pypi/pyversions/gecco-tool.svg?style=flat-square&maxAge=3600&logo=python
   :target: https://pypi.org/project/gecco-tool/#files

.. |Wheel| image:: https://img.shields.io/pypi/wheel/gecco-tool?style=flat-square&maxAge=3600&logo=data:image/png;base64,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
   :target: https://pypi.org/project/gecco-tool/#files


Overview
--------

GECCO is a fast and scalable method for identifying putative novel Biosynthetic
Gene Clusters (BGCs) in genomic and metagenomic data using Conditional Random
Fields (CRFs).

GECCO is developed in the `Zeller group <https://www.embl.de/research/units/scb/zeller/index.html>`_
and is part of the suite of computational microbiome analysis tools hosted
at `EMBL <https://www.embl.de>`_.


Quickstart
----------

Setup
^^^^^

GECCO is implemented in `Python <https://www.python.org/>`_, and supports
`all versions <https://endoflife.date/python>`_ from Python 3.6. Install
GECCO with ``pip``:

.. code:: console

  $ pip install gecco-tool

Or with Conda, using the `bioconda` channel:

.. code::

  $ conda install -c bioconda gecco


Predictions
^^^^^^^^^^^

GECCO works with DNA sequences, and loads them using Biopython, allowing it
to support a `large variety of formats <https://biopython.org/wiki/SeqIO>`_,
including the common FASTA and GenBank files.

Run a prediction on a FASTA file named ``sequence.fna`` and output the
predictions to the current directory:

.. code:: console

   $ gecco -v run --genome sequence.fna


Output
^^^^^^

GECCO will create the following files once done (using the same prefix as
the input file):

- ``{sequence}.features.tsv``: The *features* file, containing the identified
  proteins and domains in the input sequences.
- ``{sequence}.clusters.tsv``: If any were found, a *clusters* file, containing
  the coordinates of the predicted clusters, along their putative biosynthetic
  type.
- ``{sequence}_cluster_{N}.gbk``: If any were found, a GenBank file per cluster,
  containing the cluster sequence annotated with its member proteins and domains.


Reference
---------

GECCO can be cited using the following preprint:

**Accurate de novo identification of biosynthetic gene clusters with GECCO**.
Laura M Carroll, Martin Larralde, Jonas Simon Fleck, Ruby Ponnudurai, Alessio Milanese, Elisa Cappio Barazzone, Georg Zeller.
bioRxiv 2021.05.03.442509; `doi:10.1101/2021.05.03.442509 <https://doi.org/10.1101/2021.05.03.442509>`_


Feedback
--------

Contact
^^^^^^^

If you have any question about GECCO, if you run into any issue, or if you
would like to make a feature request, please create an
`issue in the GitHub repository <https://github.com/zellerlab/GECCO/issues/new>`_.
You can also directly contact `Martin Larralde via email <mailto:martin.larralde@embl.de>`_.

Contributing
^^^^^^^^^^^^

If you want to contribute to GECCO, please have a look at the
contribution guide first, and feel free to open a pull
request on the `GitHub repository <https://github.com/zellerlab/GECCO>`_.


Documentation
-------------

Guides
^^^^^^

.. toctree::
    :maxdepth: 1

    Installation <install>
    Integrations <integrations>
    Training <training>
    Contributing <contributing>

Library
^^^^^^^

.. toctree::
   :maxdepth: 1

   API reference <api/index>
   Changelog <changes>


License
-------

GECCO is released under the
`GNU General Public License v3 <https://choosealicense.com/licenses/gpl-3.0/>`_
*or later*, and is fully open-source. The ``LICENSE`` file distributed with
the software contains the complete license text.


About
-----

GECCO is developped by the `Zeller group <https://www.embl.de/research/units/scb/zeller/index.html>`_
at the European Molecular Biology Laboratory in Heidelberg. The following
individuals contributed to the development of GECCO:

- `Laura M. Carroll <https://github.com/lmc297>`_
- `Martin Larralde <https://github.com/althonos>`_
- `Jonas S. Fleck <https://github.com/joschif>`_
