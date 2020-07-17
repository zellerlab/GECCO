API Reference
==============

.. toctree::
   :hidden:

   model
   orf
   hmmer
   crf
   refine
   knn


Data Model
----------

.. currentmodule:: gecco.model

.. autosummary::
   :nosignatures:

   Hmm
   Strand
   Domain
   Protein
   Gene
   Cluster


ORF Extraction
--------------

.. currentmodule:: gecco.orf

.. autosummary::
   :nosignatures:

   ORFFinder
   PyrodigalFinder


Domain Annotation
-----------------

.. currentmodule:: gecco.hmmer

.. autosummary::
   :nosignatures:

   HMMER

BGC Detection
-------------

.. currentmodule:: gecco.crf

.. autosummary::
   :nosignatures:

   ClusterCRF


BGC Extraction
--------------

.. currentmodule:: gecco.refine

.. autosummary::
   :nosignatures:

   ClusterRefiner


Type Prediction
---------------

.. currentmodule:: gecco.knn

.. autosummary::
   :nosignatures:

   ClusterKNN
