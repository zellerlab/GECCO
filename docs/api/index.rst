API Reference
==============

.. toctree::
   :hidden:

   model
   orf
   hmmer
   interpro
   crf
   refine
   types


Data Model
----------

.. currentmodule:: gecco.model

.. autosummary::
   :nosignatures:

   ProductType
   Strand
   Domain
   Protein
   Gene
   Cluster
   ClusterTable
   FeatureTable


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

   PyHMMER

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

.. currentmodule:: gecco.types

.. autosummary::
   :nosignatures:

   TypeBinarizer
   TypeClassifier


InterPro Metadata
-----------------

.. currentmodule:: gecco.interpro

.. autosummary::
   :nosignatures:

   InterPro
   InterProEntry
