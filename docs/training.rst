Training
========

Requirements
^^^^^^^^^^^^

In order to train GECCO, you will need to have it installed with the *train*
dependencies:

.. code-block:: console

  $ git clone https://github.com/zellerlab/GECCO
  $ pip install ./GECCO[train]

This will install additional Python packages, such as `pandas <https://pandas.pydata.org/>`_
which is needed to process the feature tables, or `fisher <https://pypy.org/project/fisher>`_
which is used to select the most informative features.


Resources
^^^^^^^^^

You are free to use sequences coming from anywhere as your _positive_ regions.
GECCO was trained to detect Biosynthetic Gene Clusters, so the
`MIBiG <https://mibig.secondarymetabolites.org/>`_ database was used; you need
to get the FASTA file from MIBiG containing all BGC proteins. It can be found
`on the Download page of MIBiG <https://mibig.secondarymetabolites.org/download>`_.

You also need to have some _negative_ regions, such as bacterial genomes, and
make sure they are free of _positives_ (e.g., they must not contain any BGCs).
A common approach is to download a bunch of complete genomes from the
`ENA <https://www.ebi.ac.uk/ena/browser/home>`_ and to removes the regions where
positives are detected (e.g. removing BGCs found by `antiSMASH <https://antismash.secondarymetabolites.org/>`_).


Features
^^^^^^^^

GECCO does not train on sequences directly, but on feature tables. To obtain
a feature table from the sequences, the ``gecco annotate`` subcommand can be
used:

.. code-block:: console

  $ gecco annotate --genome reference_genome.fna -o nobgc
  $ gecco annotate --mibig mibig_prot_seqs_xx.faa -o bgc

Use the `--hmm` flag with the path to a different HMM file to use other HMMs than
the internal ones from GECCO (a subset of Pfam v33.1 and Tigrfam v15.0). **Make
sure the `--mibig` input flag is used when annotating MIBiG sequences to ensure
the metadata are properly extracted from the sequence ID of each protein in the
input file.**

_This step takes a long time, count about 5 minutes to annotate 1M amino-acids
with Pfam alone._


Embedding
^^^^^^^^^

When both the BGC and the non-BGC sequences have been annotated, merge them into
a continuous table to train the CRF on using the ``gecco embed`` subcommand:

.. code-block:: console

  $ gecco embed --bgc bgc/*.features.tsv --nobgc nobgc/*.features.tsv -o embedding.tsv


Fitting
^^^^^^^

To train the model with default parameters, you must have built a feature table
and a cluster table similar to the ones created by GECCO after a prediction.
Once it's done, use the following command to fit the weights of the CRF, and to
train the BGC type classifier:

.. code-block:: console

  $ gecco train --features embedding.features.tsv --clusters embedding.clusters.tsv -o model

Use the ``--select`` flag to select a fraction of most informative features
before training to reduce the total feature set (for instance, use ``--select 0.3``
to select the 30% features with the lowest Fisher *p-value*).

The command will output a folder (named ``model`` here) containing all the data
files necessary for predicting probabilities with ``gecco run``:

.. code-block:: console

  $ gecco run --model model ...
