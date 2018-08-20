.. PEDIA-workflow documentation master file, created by
   sphinx-quickstart on Thu May 10 16:09:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


==============================================
PEDIA-workflow usage and program documentation
==============================================

This application implements a workflow for exome variant prediction of human
genetics cases on the basis of variant, feature and photographic gestalt
information.  In order to extend the available data, a large number of cases
without exome sequencing results are used for simulation of exome simulation to
train and test the PEDIA classifier.

Gene predictions are sourced from Face2Gene and Phenomizer by translating
predicted Syndromes to genes. Detected pathogenic variants are cleaned using
Mutalyzer and converted to VCF with Jannovar.

This document contains instructions for general installation and usage of the
pipeline, but also deeper implementation details for those interested in
understanding and changing the source code.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   usage
   implementation



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
