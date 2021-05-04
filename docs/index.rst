GECCO
=====

*Biosynthetic Gene Cluster prediction with Conditional Random Fields.*


|GitLabCI| |Coverage| |License| |Source| |Issues| |Preprint|

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


Predictions
^^^^^^^^^^^

GECCO works with DNA sequences, and loads them using Biopython, allowing it
to support a `large variety of formats <https://biopython.org/wiki/SeqIO>`_,
including the common FASTA and GenBank files.

Run a prediction on a FASTA file named ``sequence.fna`` and output the
predictions to the current directory:

.. code:: console

   $ gecco run --genome sequence.fna


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
