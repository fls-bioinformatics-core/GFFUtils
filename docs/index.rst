GFFUtils: utilities for handling GFF and GTF files
==================================================

Overview
********

The ``GFFUtils`` package provides a small set of utility programs
for working with GFF and GTF files, specifically:

 * ``gff_cleaner``: perform "cleaning" operations on a GFF file
 * ``gff_annotation_extractor``: combine and annotate feature
   counts (e.g. from ``htseq-count``) with data from a GFF file
 * ``gtf_extract``: extract selected data items from a GTF file
 * ``gtf2bed``: convert GTF file to BED format

.. warning::

   The old names for the utilities (``GFFcleaner``,
   ``GFF3_Annotation_Extractor`` and ``GTF_extract``) are
   still supported, but will be removed in a future release.

Installation
************

To install the latest version of ``GFFUtils``, download the latest
release from as a ``tar.gz`` file from:

https://github.com/fls-bioinformatics-core/GFFUtils/releases

For example for version 0.10.3, download ``GFFUtils-0.10.3.tar.gz``.

Unpack the code using:

::

   tar xzf GFFUtils-0.10.3.tar.gz

which will unpack into a new directory called e.g.
``GFFUtils-0.10.3``.

It is recommended to install the code into a Python virtual
environment, which you can create by doing:

::

   virtualenv venv
   source venv/bin/activate
   pip install -r ./GFFUtils-0.10.3/requirements.txt
   pip install ./GFFUtils-0.10.3/

To install the developmental code directly from GitHub:

::

    pip install -r https://raw.githubusercontent.com/fls-bioinformatics-core/GFFUtils/devel/requirements.txt
    pip install git+https://github.com/fls-bioinformatics-core/GFFUtils.git@devel

.. note::

   ``GFFUtils`` is currently supported under the following Python
   versions:

   * 2.7
   * 3.5
   * 3.6
   * 3.7
   * 3.8

   but support for Python 2.7 is likely to be dropped in the near
   future.

Contents
********

.. toctree::
   :maxdepth: 2

   gff_cleaner
   gff_annotation_extractor
   gtf_extract
   gtf2bed
   extras

Version history
***************

See the :doc:`CHANGELOG <changes>`.

Additional information
**********************

See http://www.sanger.ac.uk/resources/software/gff/spec.html for
details of the GFF format.



Credits
*******

These utilities have been developed by Peter Briggs with input from
Leo Zeef, to support the activities of the Bioinformatics Core Facility
(BCF) in the Faculty of Biology Medicine and Health (FBMH) at the
University of Manchester (UoM).

