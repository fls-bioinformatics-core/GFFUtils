``GFFcleaner``: clean up GFF files
==================================

Overview
--------

The ``GFFcleaner`` utility performs various manipulations on a GFF
file to "clean" it.

Usage and options
-----------------

General usage syntax::

     GFFcleaner [OPTIONS] <file>.gff

Options:

.. cmdoption:: --version

   show program's version number and exit

.. cmdoption:: -h, --help

   show the help message and exit

.. cmdoption:: -o OUTPUT_GFF

   Name of output GFF file (default is ``<file>_clean.gff``)

.. cmdoption:: --prepend=PREPEND_STR

   String to prepend to seqname in first column

.. cmdoption:: --clean

   Perform all the 'cleaning' manipulations on the input
   data (equivalent to specifying all of ``--clean-score``,
   ``--clean-replace-attributes``, ``--clean-exclude-attributes``
   and ``--clean-group-sgds``)

.. cmdoption:: --clean-score

   Replace ``Anc_*`` and blanks in the ``score`` field with
   zeroes

.. cmdoption:: --clean-replace-attributes

   Replace ``ID``, ``Gene``, ``Parent`` and ``Name`` attributes
   with the value of the ``SGD`` attribute, if present

.. cmdoption:: --clean-exclude-attributes

   Remove the ``kaks``, ``kaks2`` and ``ncbi`` attributes (to
   remove arbitrary attributes, see the ``--remove-attribute=...``
   option)

.. cmdoption:: --clean-group-sgds

   Group features with the same ``SGD`` by adding unique
   numbers to the ``ID`` attributes; ``ID``s will have the form
   ``CDS:<SGD>:<n>`` (where ``n`` is a unique number for a
   given SGD)

.. cmdoption:: --report-duplicates

   Report duplicate ``SGD`` names and write list to
   ``<file>_duplicates.gff`` with line numbers, chromosome,
   start coordinate and strand.

.. cmdoption:: --resolve-duplicates=MAPPING_FILE

   Resolve duplicate ``SGD``s by matching against 'best'
   genes in the supplied mapping file; other non-matching genes
   are discarded and written to ``<file>_discarded.gff``.

.. cmdoption:: --discard-unresolved

   Discard any unresolved duplicates, which are written
   to ``<file>_unresolved.gff``.

.. cmdoption:: --insert-missing=GENE_FILE

   Insert genes from gene file with ``SGD`` names that don't
   appear in the input GFF. If ``GENE_FILE`` is blank ('='s
   must still be present)  then the mapping file supplied
   with the ``--resolve-duplicates`` option will be used
   instead.

.. cmdoption:: --add-exon-ids

   For exon features without an ``ID`` attribute, construct
   and insert an ID of the form ``exon:<Parent>:<n>``
   (where ``n`` is a unique number).

.. cmdoption:: --add-missing-ids

   For features without an ``ID`` attribute, construct and
   insert a generated ID of the form ``<feature>:<Parent>:<n>``
   (where ``n`` is a unique number).

.. cmdoption:: --no-percent-encoding

   Convert encoded attributes to the correct characters
   in the output GFF.

   .. warning::

   This may result in a non-cannonical GFF that can't be read
   correctly by this or other programs.

.. cmdoption:: --remove-attribute=RM_ATTR

   Remove attribute ``RM_ATTR`` from the list of attributes
   for all records in the GFF file (can be specified
   multiple times)

.. cmdoption:: --strict-attributes

   Remove attributes that don't conform to the ``KEY=VALUE``
   format

.. cmdoption:: --debug

   Print debugging information

Output files
------------

 * ``<file>_clean.gff``: 'cleaned' version of input
 * ``<file>_duplicates.txt``: list of duplicated ``SGD`` names
   and the lines they appear on in the input file, along with
   chromosome, start coordinate and strand
 * ``<file>_discarded.gff``: genes rejected by
   ``--resolve-duplicates``
 * ``<file>_unresolved.gff``: unresolved duplicates rejected by
   ``--discard-unresolved``

Usage recipe
------------

The following steps outline the procedure for using the program,
with each step being run on the output from the previous one:

1. **Clean the chromosome names in the file by adding a prefix
   (`--prepend` option)**

   Creates a copy of the input file with the chromosome names
   updated with a specified prefix.

   E.g. ``--prepend=chr`` will add ``chr`` to the start of each
   chromosome name in the file, which is useful if the chromosome
   is denoted by a number and needs the prefix for consistency with
   a mapping file.

2. **Clean the GFF score and attribute data (`--clean` options)**

   The "clean" options perform the following operations:

   * ``--clean-score``: the data in the score column is cleaned up
     by replacing ``Anc_*`` and blanks with '0's.

   The attribute field of the GFF can contain various
   semicolon-separated key-value pairs:

   * ``--clean-replace-attributes``: if one of these is a non-blank
     ``SGD`` then the ``Gene``, ``Parent`` and ``Name`` values are
     updated to be the same as the ``SGD`` name.
   * ``--clean-exclude-attributes``: attributes called ``kaks``,
     ``kaks2`` and ``ncbi`` are removed (n.b. to remove arbitrary
     attributes, use the more general ``--remove-attribute=...``
     option).

   If multiple features share the same SGD name then
   ``--clean-replace-attributes`` can result in them also sharing the
   same ID; to deal with this:

   * ``--clean-group-sgds``: update the ID attribute to group
     neighbouring lines that have the same ``SGD`` (see
     :ref:`sgd_grouping` below).

     A single `--clean` can be specified which performs all these operations
     automatically.

3. **Detect duplicate SGDs (`--report-duplicates` option)**

   Report duplicate SGD names found in the input file.

   This option writes a list of the duplicates to a 'duplicates' file.

   It also reports the number of 'trivial' duplicates, i.e. lines having
   the same ``SGD`` because they are part of the same gene.

4. **Resolve duplicate SGDs using a mapping file (`--resolve-duplicates`
   option)**

   Attempt to resolve duplicates by referring to a list of "best" genes
   given in a mapping file. For each duplicated name the resolution
   procedure is:

   * Find mapping gene(s) with the same name
   * For each mapping gene, keep duplicates which match chromosome, strand
     and which overlap with the start and end of the gene (see
     :ref:`overlap_criteria` below). For ``SGD`` groups the mapping gene must
     overlap the whole group for it to match; mapping genes and duplicates
     which don't have matches are removed from the process.
   * At the end of the matching procedure the duplication is resolved if
     there is one ``SGD`` (or ``SGD`` group) matched to one mapping gene.
     Otherwise the duplication remains unresolved.

   When duplicates are resolved, the non-matching duplicates are discarded;
   otherwise by default all unresolved duplicates are kept. However if the
   ``--discard-unresolved`` option is also specified then all unresolved
   duplicates are removed before output; the ``--insert-missing`` option
   can then be used to add them back in.

   Note that the ``--discard-unresolved`` option cannot get rid of 'trivial'
   duplicates (i.e. lines having the same SGD because they are part of the
   same gene).

5. **Add missing genes (`--insert-missing` option)**

   Adds genes from a list of "best" genes given in a mapping file which
   have names not found in the input GFF.

.. _`sgd_grouping`:

SGD grouping
------------

As part of setting the ``ID`` attribute of GFF lines, the "clean" option
also attempts to group neighbouring lines which have the same ``SGD`` name.

The ``ID`` attribute is updated to the form::

    ID=CDS:<sgd_name>:<i>

where ``<sgd_name>`` is a gene or transcript name (e.g. ``YEL0W``) and
``<i>`` is an integer index which starts from 1. Groupings are indicated
by subsequent lines having the same ``<sgd_name>`` but monotonically
increasing indices, for example::

    chr1   Test   CDS   34525   35262   0   -   0   ID=CDS:YEL0W:1;SGD=YEL0W
    chr1   Test   CDS   35823   37004   0   -   0   ID=CDS:YEL0W:2;SGD=YEL0W
    chr1   Test   CDS   38050   38120   0   -   0   ID=CDS:YEL0W:3;SGD=YEL0W
    chr1   Test   CDS   39195   39569   0   -   0   ID=CDS:YEL0W:4;SGD=YEL0W

When determining a grouping the program looks ahead from each line for
subsequent lines (up to five) which have the same SGD value. So groupings
can also accommodate "breaks", for example::

    chr1   Test   CDS   34525   35262   0   -   0   ID=CDS:YEL0W:1;SGD=YEL0W
    chr1   Test   CDS   35823   37004   0   -   0   ID=CDS:YEL0W:2;SGD=YEL0W
    chr1   Test   CDS   38050   38120   0   -   0   ID=CDS:YEL0X:1;SGD=YEL0X
    chr1   Test   CDS   39195   39569   0   -   0   ID=CDS:YEL0W:3;SGD=YEL0W

Mapping file format
-------------------

The mapping file is a tab-delimited text file with lines of the form::

    name   chr   start   end    strand

``<name>`` is used to match against the ``SGD`` names in the input GFF file.

.. _`overlap_criteria`:

Overlap criteria
----------------

Aside from matching chromosome and strand, one of the criteria for a
mapping gene to match a duplicate from the GFF file is that the two must
overlap.

An overlap is counted as the duplicate from the GFF having start/end
positions such that it lies inside the start/end positions of the mapping
gene extended by 1kb i.e. between ``start - 1000`` and ``end + 1000``.
