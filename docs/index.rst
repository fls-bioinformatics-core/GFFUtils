GFFUtils: utilities for handling GFF and GTF files
==================================================

Overview
********

The ``GFFUtils`` package provides a small set of utility programs
for working with GFF and GTF files, specifically:

 * ``GFFcleaner``: perform "cleaning" operations on a GFF file
 * ``GFF3_Annotation_Extractor``: combine and annotate feature
   counts (e.g. from ``htseq-count``) with data from a GFF file
 * ``GTF_extract``: extract selected data items from a GTF file
 * ``gtf2bed``: convert GTF file to BED format

Installation
************

To install from github, do::

    pip install -r https://raw.githubusercontent.com/fls-bioinformatics-core/GFFUtils/refactor/requirements.txt
    pip install git+https://github.com/fls-bioinformatics-core/GFFUtils.git


Contents
********

.. toctree::
   :maxdepth: 2

   install
   GFFcleaner
   GFF3_Annotation_Extractor
   GTF_extract
   gtf2bed
   extras

Additional information
**********************

See http://www.sanger.ac.uk/resources/software/gff/spec.html for
details of the GFF format.

Credits
*******

These utilities have been developed by Peter Briggs with input from
Leo Zeef, to support the activities of the Bioinformatics Core Facility
(BCF) in the Faculty of Life Sciences (FLS) at the University of
Manchester (UoM).


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

