GFFUtils
========

Utilities for working with GFF and GTF files:

* ``gff_cleaner.py``: perform various "cleaning" manipulations on a GFF file
* ``gff_annotation_extractor``: combine and annotate feature counts with data
  from a GFF or GTF file
* ``gft_extract``: extract selected data items from a GTF file
* ``gtf2bed``: convert GTF contents to BED format

Full documentation is available at http://gffutils.readthedocs.org/

Installation
------------

To install the latest version of ``GFFUtils``, download the latest
release from as a ``tar.gz`` file from:

https://github.com/fls-bioinformatics-core/GFFUtils/releases

For example for version 0.10.3, download ``GFFUtils-0.10.3.tar.gz``.

Unpack the code using:

::

   tar xzf GFFUtils-0.10.3.tar.gz

which will unpack into a new directory called e.g. ``GFFUtils-0.10.3``.

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

Note that ``GFFUtils`` is currently supported for the following Python
versions:

* 2.7
* 3.5
* 3.6
* 3.7
* 3.8

but support for Python 2.7 is likely to be dropped in the near future.

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

.. image:: https://travis-ci.org/fls-bioinformatics-core/GFFUtils.png?branch=master
   :alt: Current status of TravisCI build for master branch
   :target: https://travis-ci.org/fls-bioinformatics-core/GFFUtils/builds

Credits
-------

These utilities have been developed by Peter Briggs with input from
Leo Zeef, to support the activities of the Bioinformatics Core Facility
(BCF) in the Faculty of Biology Medicine and Health (FBMH) at the
University of Manchester (UoM).
