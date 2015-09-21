GFFUtils
========

Utilities for working with GFF and GTF files:

* ``GFFcleaner.py``: perform various "cleaning" manipulations on a GFF file
* ``GFF3_Annotation_Extractor``: combine and annotate feature counts with data
   from a GFF or GTF file
* ``GTF_extract``: extract selected data items from a GTF file

Installation
------------

It is recommended to use::

    pip install -r requirements.txt
    pip install .

from within the top-level source directory to install the package.

To use the package without installing it first you will need to add the
directory to your ``PYTHONPATH`` environment.

To install directly from github using ``pip``::

    pip install git+https://github.com/fls-bioinformatics-core/GFFUtils.git

(In either of these latter two cases you will also need to install the
``genomics-bcftbx`` package from
https://github.com/fls-bioinformatics-core/genomics)

Documentation
-------------

Documentation based on ``sphinx`` is available under the ``docs`` directory.

To build do either::

    python setup.py sphinx_build

or::

    cd docs
    make html

both of which create the documentation in the ``docs/build`` subdirectory.

Running Tests
-------------

The Python unit tests can be run using::

    python setup.py test

Note that this requires the ``nose`` package.

In addition the tests are run via TravisCI whenever this GitHub repository
is updated:

.. image:: https://travis-ci.org/fls-bioinformatics-core/GFFUtils.png?branch=mas
ter
   :alt: Current status of TravisCI build for master branch
   :target: https://travis-ci.org/fls-bioinformatics-core/genomics/builds

Credits
-------

These utilities have been developed by Peter Briggs with input from
Leo Zeef, to support the activities of the Bioinformatics Core Facility
(BCF) in the Faculty of Life Sciences (FLS) at the University of
Manchester (UoM).
