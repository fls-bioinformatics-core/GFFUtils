``GFF3_Annotation_Extractor``: annotate gene feature data
=========================================================

Overview
--------

``GFF3_Annotation_Extractor`` takes gene feature data (for example the
output from one or more runs of the HTSeq-count program) and combines
it with data about each feature's parent gene from a GFF file.

By default the program takes a single tab-delimited input file where the
first column contains feature IDs, and appends data about the feature's
parent gene.

In 'htseq-count' mode, one or more ``htseq-count`` output files should
be provided as input, and the program will write out the data about the
feature's parent gene appended with the counts from each input file.

To generate the feature count files using ``htseq-count`` do e.g.::

    htseq-count --type=exon -i Parent <file>.gff <file>.sam

which returns counts of each exon against the name of that exon's parent.

``GFF3_Annotation_Extractor`` will match up the exon parent with its
parent gene and output the counts against gene names.

.. note::

   ``GFF3_Annotation_Extractor`` can also be used with GTF input.

Usage and options
-----------------

General usage syntax::

    GFF3_Annotation_Extractor.py OPTIONS <file>.gff FEATURE_DATA

or::

    GFF3_Annotation_Extractor.py --htseq-count OPTIONS <file>.gff FEATURE_COUNTS [FEATURE_COUNTS2 ...]

if working with ``htseq-count`` data.

Options:

.. cmdoption:: --version

   show program's version number and exit

.. cmdoption:: -h, --help

   show the help message and exit

.. cmdoption:: -o OUT_FILE

   specify output file name

.. cmdoption:: -t FEATURE_TYPE, --type=FEATURE_TYPE

   feature type listed in input count files (default
   ``exon``; if used in conjunction with ``--htseq-count``
   option then should be the same as that specified when
   running htseq-count)

.. cmdoption:: -i ID_ATTRIBUTE, --id-attribute=ID_ATTRIBUTE

   explicitly specify the name of the attribute to get
   the feature IDs from (defaults to ``ID`` for GFF input,
   ``gene_id`` for GTF input)

.. cmdoption:: --htseq-count

   htseq-count mode: input is one or more output
   ``FEATURE_COUNT`` files from the ``htseq-count`` program

Output files
------------

* ``<basename>_annot.txt``: the feature data annotated with data
  for each parent gene.
* ``<basename>_annot_stats.txt``: the counts of "ambiguous",
  "two_low_aQual" etc from each log (htseq-count mode only).

