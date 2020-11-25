Version History and Changes
===========================

---------------------------
Version 0.12.0 (2020-11-25)
---------------------------

* ``gff_annotation_extractor``: fix bug with ``-t``/``--type`` option
  which previous didn't do anything
  (`PR #49 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/49>`_)
* ``gff_annotation_extractor``: fix bug when annotating non-htseq-count
  feature files which didn't have a header
  (`PR #48 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/48>`_)
* ``gff_annotation_extractor``: update the documentation to clarify
  function and usage
  (`PR #47 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/47>`_)
  

---------------------------
Version 0.11.0 (2020-09-02)
---------------------------

* Updated licence to Academic Free Licence 3
* Standardise the utility names to ``gff_annotation_extractor``,
  ``gff_cleaner`` and ``gtf_extract`` (NB the old names are still
  available but will be dropped in a future release)
  (`PR #41 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/41>`_)
* Updated to support Python 3 (versions 3.5, 3.6, 3.7 and 3.8);
  Python 2.7 is still supported but will be dropped in a future
  release
  (`PR #34 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/34>`_)
* ``gtf2bed``: new option ``-o`` allows output file to be specified
  (`PR #28 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/28>`_)

---------------------------
Version 0.10.3 (2019-09-10)
---------------------------

Bugfix patch release.

* ``GFFcleaner:`` fix bug with ``--clean`` option when GFF attributes
  don't include SGD values (`PR #12 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/12>`_)

---------------------------
Version 0.10.2 (2019-09-09)
---------------------------

Bugfix patch release.

* ``GFFcleaner:``: fix bug with ``--clean-score`` option when GFF has
  score values which are not blank, don't start with ``Anc_``, or are
  non-zero (`PR #11 <https://github.com/fls-bioinformatics-core/GFFUtils/pull/11>`_)
  
---------------------------
Version 0.10.1 (2019-02-11)
---------------------------

Patch release to update installation documentation (no other functional
changes or bug fixes)

---------------------------
Version 0.10.0 (2018-10-04)
---------------------------

General:

* Documentation updated to indicate that only Python 2 is supported;
  ``setup.py`` also includes python_requires keyword which is valid
  for ``pip`` >= 9.

``GTF_extract``:

* Fix the broken ``-o`` option.
* New option ``-k``/``--keep-header`` will copy the header of the input
  file to the output.
