``GTF_extract``: extract data items from GTF/GFF
================================================

Overview
--------

The ``GTF_extract`` utility extracts selected data items from a
GTF file and output in tab-delimited format.

.. note::

   The program can also operate on GFF files provided the ``--gff``
   option is specified.

Usage and options
-----------------

General usage syntax::

    GTF_extract.py OPTIONS <gft_file>

Options:

.. cmdoption:: --version

   show program's version number and exit

.. cmdoption:: -h, --help

   show the help message and exit

.. cmdoption:: -f FEATURE_TYPE, --feature=FEATURE_TYPE

   only extract data for lines where feature is FEATURE_TYPE

.. cmdoption:: --fields=FIELD_LIST

   comma-separated list of fields to output in tab-delimited format
   for each line in the GTF, e.g. ``chrom,start,end``.

   Fields can either be a GTF field name (i.e. ``chrom``, ``source``,
   ``feature``, ``start``, ``end``, ``score``, ``strand`` and
   ``frame``), or the name of an attribute (e.g. ``gene_name``,
   ``gene_id`` etc).

   Data items are output in the order they appear in ``FIELD_LIST``.
   If a field doesn't exist for a line then ``'.'`` will be output as
   the value.

.. cmdoption:: -o OUTFILE

   write output to OUTFILE (default is to write to stdout)

.. cmdoption:: --gff

   specify that the input file is GFF rather than GTF format

Output
------

The program outputs a tab-delimited line of data for each matching line
found in the input GTF file; the data items in the line are those
specified by the ``--fields`` option (or else all data items, if no fields
were specified).

For example, for ``--fields=chrom,start,end,strand``, the GTF line::

    chr1	HAVANA	gene	11869	14412	.	+	.	gene_id "ENSG00000223972.4" ...

will produce the output::

    chr1	11869	14412	+

By default the output of the program is written to stdout; use the
``-o`` option to direct the output to a named file instead.
