``gff_annotation_extractor``: annotate gene feature data
========================================================

Overview
--------

``gff_annotation_extractor`` takes gene feature data (for example the
output from one or more runs of the HTSeq-count program) and combines
it with data about each feature's parent gene from a GFF file.

By default the program takes feature data from a single tab-delimited
input file where the first column contains feature IDs, and outputs an
updated copy of the file with data about the feature's parent feature
and parent gene appended to each line.

In 'htseq-count' mode, one or more ``htseq-count`` output files should
be provided as input; the program will write out the data about the
feature's parent feature and parent gene appended with the counts from
each input file.

By default feature IDs from the feature data files are matched to
the first record in the input GFF where the ``ID`` attribute of that
record is the same (a different attribute can be specified using the
``-i`` option). All records are considered regardless of the feature
type, unless the ``-t`` option is used to restrict the records to
just those with the specified feature type (this may be required in
'htseq-count' mode).

The parent gene is located by recursively looking up records where
the ``ID`` attribute matches the ``Parent`` attribute, until a
gene record is found.

.. note::

   ``gff_annotation_extractor`` can also be used with GTF input,
   in which case the feature IDs are matched using the ``gene_id``
   attribute by default. Only ``gene`` feature types are considered
   when using GTF data.

Usage and options
-----------------

General usage syntax:

::

    gff_annotation_extractor OPTIONS <file>.gff FEATURE_DATA

Usage in 'htseq-count' mode:

::

    gff_annotation_extractor --htseq-count OPTIONS <file>.gff FEATURE_COUNTS [FEATURE_COUNTS2 ...]

Options:

.. cmdoption:: --version

   show program's version number and exit

.. cmdoption:: -h, --help

   show the help message and exit

.. cmdoption:: -o OUT_FILE

   specify output file name

.. cmdoption:: -t FEATURE_TYPE, --type=FEATURE_TYPE

   restrict feature records to this type when matching
   features from input count files; if used in conjunction
   with ``--htseq-count`` then should be the same as that
   specified when running htseq-count (default: include all
   feature records)

.. cmdoption:: -i ID_ATTRIBUTE, --id-attribute=ID_ATTRIBUTE

   explicitly specify the name of the attribute to get
   the feature IDs from (defaults to ``ID`` for GFF input,
   ``gene_id`` for GTF input)

.. cmdoption:: --htseq-count

   htseq-count mode: input is one or more output
   ``FEATURE_COUNT`` files from the ``htseq-count`` program

'htseq-count' mode
------------------

To generate the feature count files using ``htseq-count`` do e.g.:

::

    htseq-count --type=exon -i Parent <file>.gff <file>.sam

which returns counts of each exon against the name of that exon's
parent.

``gff_annotation_extractor`` should then be run using the same
value for the ``--type`` option:

::

   gff_annotation_extractor --htseq-count --type=exon <file>.gff <counts>.out

Output files
------------

``gff_annotation_extractor`` always produces a copy of the feature
data annotated with data for each parent gene. By default this will
be called ``<basename>_annot.txt``; use the ``-o`` option to specify
a different name.

The annotation consists of the following fields:

* ``exon_parent``: ID for the parent feature
* ``feature_type_exon_parent``: type for the parent feature
* ``gene_ID``: ID for the gene the feature belongs to
* ``gene_name``: name of the gene (from the ``Name`` attribute for
  GFF, or ``gene_name`` attribute for GTF)
* ``chr``: chromosome of the gene
* ``start``: start position of the gene
* ``end``: end position of the gene
* ``strand``: strand for the gene
* ``gene_length``: gene length
* ``locus``: string consisting of ``<chr>:<start>-<end>``
* ``description``: text from the gene's ``description`` attribute

In the default mode these fields are appended to each line from
the input feature file; in 'htseq-count' mode each line in the
annotation file consists of these fields, with the counts from
each ``htseq-count`` file appended.

If a parent gene cannot be located for a feature then the
annotation for that feature will be empty.

In 'htseq-count' mode an additonal file called
``<basename>_annot_stats.txt`` is also produced with the counts
of "ambiguous", "two_low_aQual" etc from each log.

Warnings and errors
-------------------

The following is a non-exhaustive list of the warnings and errors
that ``gff_annotation_extractor`` can produce, along with a brief
description and possible cause:

* ``Unable to locate parent data for feature '...'``: indicates IDs
  in the feature files for which no matching records can be located
  in the input GFF. In this case the output annotation will be blank.
  Check that the input feature file consists of tab-delimited data.

* ``Multiple parents found on line ...``: indicates that a record
  matching a feature ID has a ``Parent`` attribute which contains
  multiple comma-separated IDs. In this case it may not be possible
  to locate the parent gene for the feature.

* ``No identifier attribute (...) on line ...``: indicates a record
  from the input GFF with no ``ID`` attribute (or custom attribute
  supplied via ``-i`` option).

* ``No '...' attribute found on line ...``: indicates a record
  from the input GTF with no ``gene_id`` attribute (or custom
  attribute supplied via ``-i`` option).
