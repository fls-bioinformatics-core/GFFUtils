Version History and Changes
===========================

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
