Training
========


Software Requirements
---------------------

In order to train GECCO, you need to have installed with the *train* optional dependencies.
This can be done with ``pip``:

.. code-block:: console

  $ pip install gecco-tool[train]

This will install additional Python packages, such as `pandas <https://pandas.pydata.org/>`_
which is needed to process the feature tables, or `fisher <https://pypy.org/project/fisher>`_
which is used to select the most informative domains.


Domain database
---------------

GECCO needs HMM domains to use as features. Installing the ``gecco-tool`` package
will also install a subset of the Pfam database that can be used for making the
predictions. However, this subset should not be used for training, since a
different subset of domains may be selected with different training data.

You can get the latest version of Pfam (35.0 in December 2021) from the EMBL
FTP server:

.. code::

    $ wget "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz" -O Pfam35.hmm.gz

You are also free to get additional HMMs from other databases, such as
`TIGRFAM <https://www.jcvi.org/research/tigrfams>`_ or `PANTHER <http://www.pantherdb.org/panther/;jsessionid=D7BFDD605F98EC1159A5E0E77536FD76>`_,
or even to build your own HMMs, as long as they are in `HMMER <http://hmmer.org/>`_ format.


Training sequences
------------------

Regions in their genomic context
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest case for training GECCO is when you have entire genomes with regions
of interest. In that case, you can directly use these sequences, and you will
only have to prepare a cluster table with the coordinates of each positive region.


Regions from separate datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the event you don't have the genomic context available for your regions of
interest, you will have to provide a "fake" context by embedding the positive
regions into contigs that don't contain any positive.

GECCO was trained to detect Biosynthetic Gene Clusters, so the
`MIBiG <https://mibig.secondarymetabolites.org/>`_ database was used to get
the positives, with some additional filtering to remove redundancy and entries
with invalid annotations. For the negative regions, we used representative
genomes from the `proGenomes <https://progenomes.embl.de/>`_ database, and masked
known BGCs using `antiSMASH <https://antismash.secondarymetabolites.org/>`_.



Gene and Feature tables
-----------------------

GECCO does not train on sequences directly, but on gene and feature tables. 
You can build the tables yourself (see below for the expected format), but the 
easiest way to obtain a gene table and a feature table from any sequence is the 
``gecco annotate`` subcommand. To build a table from a collection of nucleotide 
sequences in ``sequences.fna`` and HMMs in ``Pfam35.hmm.gz``, use:

.. code-block:: console

  $ gecco annotate --genome sequences.fna --hmm Pfam35.hmm.gz -o .

This will output two files in the local folder, ``sequences.genes.tsv`` and 
``sequences.features.tsv`` containing the gene table and the feature table for 
the input sequences.

.. hint::

    If this step takes too long, you can also split the file containing your
    input sequences, process them independently in parallel, and combine the
    result.

The gene table is a TSV file that looks like this, with one row per gene, per
sequence:

============  ============== ===== ==== ======
sequence_id   protein_id     start end  strand
============  ============== ===== ==== ======
AFPU01000001  AFPU01000001_1     3 2555  ``+``
AFPU01000001  AFPU01000001_2  2610 4067  ``-``
         ...             ...   ...  ...    ...
============  ============== ===== ==== ======


The feature table is a TSV file that looks like this, with one row per domain,
per protein, per sequence:

============  ============== ===== ==== ====== ======= ====== ======== =========== ============ ==========
sequence_id   protein_id     start end  strand domain  hmm    i_evalue pvalue      domain_start domain_end
============  ============== ===== ==== ====== ======= ====== ======== =========== ============ ==========
AFPU01000001  AFPU01000001_1     3 2555  ``+`` PF01573 Pfam35 95.95209 0.00500                2         27
AFPU01000001  AFPU01000001_2  2610 4067  ``-`` PF17032 Pfam35  0.75971 3.961e-05             83        142
AFPU01000001  AFPU01000001_2  2610 4067  ``-`` PF13719 Pfam35  4.89304 0.000255              85         98
         ...             ...   ...  ...    ...     ...    ...      ...         ...          ...        ...
============  ============== ===== ==== ====== ======= ====== ======== =========== ============ ==========


Cluster tables
--------------

The cluster table is used to pass additional information to GECCO: the location of
each positive region in the input data, and the type of each region (if it makes
sense). You need to build this table manually, but it should be quite straightforward.

.. hint::

    If a region has more than one type, use `;` to separate the two types
    in the type column. For instance, a Polyketide/NRP hybrid cluster can be
    marked with the type ``Polyketide;NRP``.

The cluster table is a TSV file that looks like this, with one row per region:

============ ============ ====== ====== ==========
sequence_id  cluster_id   start  end    type
============ ============ ====== ====== ==========
AFPU01000001   BGC0000001 806243 865563 Polyketide
MTON01000024   BGC0001910 129748 142173    Terpene
         ...          ...    ...    ...        ...
============ ============ ====== ====== ==========

.. hint::

    If the concept of "type" makes no sense for the regions you are trying to
    detect, you can omit the ``type`` column entirely. This will effectively
    mark all the regions from the training sequences as "Unknown". Later on,
    GECCO will skip training the cluster type classifier if it detects that 
    all the clusters are of the same type.


Fitting the model
-----------------

Now that you have everything needed, it's time to train GECCO! Use the
following method to fit the CRF model and the type classifier:

.. code-block:: console

  $ gecco -vv train --genes genes.tsv --features features.tsv --clusters clusters.tsv -o model

GECCO will create a directory named ``model`` containing all the required files
to make predictions later on. 


L1/L2 regularisation
^^^^^^^^^^^^^^^^^^^^

Use the ``--c1`` and ``--c2`` flags to control the weight for the L1 and L2
regularisation, respectively. The command line defaults to *0.15* and *0.15*;
however, for training GECCO, we disabled L2 regularisation and selected
a value of *0.4* for :math:`C_1` by optimizing on an external validation dataset.


Feature selection
^^^^^^^^^^^^^^^^^

GECCO supports selecting the most informative features from the training dataset
using a simple contingency testing for the presence/absence of each domain in
the regions of interest. Reducing the number of features helps the CRF model to
get better accuracy. It also greatly reduces the time needed to make predictions
by skipping the HMM annotation step for useless domains.

Use the ``--select`` flag to select a fraction of most informative features
before training to reduce the total feature set (for instance, use ``--select 0.3``
to select the 30% features with the lowest Fisher *p-value*).

.. code-block:: console

  $ gecco train --features features.tsv --clusters clusters.tsv -o model --select 0.3

.. hint::

    You will get a warning in case you select a *p-value* threshold that is still
    too high, resulting in non-informative domains to be included in the selected
    features.


Predicting with the new model
-----------------------------

To make predictions with a model different from the one embedded in GECCO, you
will need the folder from a previous ``gecco train`` run, as well as the HMMs
used to build the feature tables in the first place.

.. code-block:: console

  $ gecco run --model model --hmm Pfam35.hmm.gz --genome genome.fa -o ./predictions/

Congratulations, you trained GECCO with your own dataset, and successfully
used it to make predictions!
