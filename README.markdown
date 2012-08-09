GFFUtils
========

GFFUtils package provides the following utilities for working with GFF files:

 * `GFFcleaner.py`: performs various "cleaning" manipulations on a GFF file
 * `GFF3_Annotation_Extractor.py`: combine and annotate feature counts with data
   from a GFF file

GFFcleaner.py
-------------

The GFFcleaner can perform various manipulations on a GFF file to "clean" it.

### Usage ###

     GFFcleaner.py [OPTIONS] <file>.gff

### Options ###

    --version             show program's version number and exit

    -h, --help            show this help message and exit

    -o OUTPUT_GFF         Name of output GFF file (default is
                          '<file>_clean.gff')

    --prepend=PREPEND_STR
                          String to prepend to seqname in first column

    --clean               Perform all the 'cleaning' manipulations on the input
                          data (equivalent to specifying all of --clean-score,
                          --clean-replace-attributes, --clean-exclude-attributes
                          and --clean-group-sgds)

    --clean-score         Replace 'Anc_*' and blanks in 'score' field with
                          zeroes

    --clean-replace-attributes
                         Replace 'ID', 'Gene', 'Parent' and 'Name' attributes
                         with the value of the SGD attribute, if present

    --clean-exclude-attributes
                          Remove the 'kaks', 'kaks2' and 'ncbi' attributes

    --clean-group-sgds    Group features with the same SGD by adding unique
                          numbers to the 'ID' attributes; IDs will have the form
                          'CDS:<SGD>:<n>' (where n is a unique number for a
                          given SGD)

    --report-duplicates   Report duplicate SGD names and write list to
                          <file>_duplicates.gff with line numbers, chromosome,
                          start coordinate and strand.

    --resolve-duplicates=MAPPING_FILE
                          Resolve duplicate SGDs by matching against 'best'
                          genes in the supplied mapping file; other non-matching
                          genes are discarded and written to
                          <file>_discarded.gff.

    --discard-unresolved  Discard any unresolved duplicates, which are written
                          to <file>_unresolved.gff.

    --insert-missing=GENE_FILE
                          Insert genes from gene file with SGD names that don't
                          appear in the input GFF. If GENE_FILE is blank ('='s
                          must still be present) then the mapping file supplied
                          with the --resolve-duplicates option will be used
                          instead.

    --add-exon-ids        For exon features without an ID attribute, construct
                          and insert an ID of the form 'exon:<Parent>:<n>'
                          (where n is a unique number).

    --add-missing-ids     For features without an ID attribute, construct and
                          insert a generated ID of the form
                          '<feature>:<Parent>:<n>' (where n is a unique number).

    --no-percent-encoding
                          Convert encoded attributes to the correct characters
                          in the output GFF. WARNING this may result in a non-
                          cannonical GFF that can't be read correctly by this or
                          other programs.

    --debug               Print debugging information

    --test                Run unit tests

### Output files ###

 *  `<file>_clean.gff` 'cleaned' version of input

 *  `<file>_duplicates.txt` list of duplicated SGD names and the lines
     they appear on in the input file, chromosome, start coordinate and strand

 *  `<file>_discarded.gff` genes rejected by `--resolve-duplicates`

See <http://www.sanger.ac.uk/resources/software/gff/spec.html> for details of the GFF
format.

### run_cleanup.sh ###

This is a shell script which automatically runs all the steps outlined in the
_Usage recipe_.

Usage: `run_cleanup.sh <gff_file> [<mapping_file>]`

Note that duplicate resolution and missing gene insertion cannot be performed without
a mapping file.

### Usage recipe ###

The following steps outline the procedure for using the program, with each step
being run on the output from the previous one:

1.   **Clean the chromosome names in the file by adding a prefix (`--prepend option`)**

     Creates a copy of the input file with the chromosome names updated with a specified
     prefix.

     E.g. `--prepend=chr` will add `chr` to the start of each chromosome name in the
     file. 

     This is useful if the chromosome is denoted by a number and needs the prefix for
     consistency with a mapping file.

2.   **Clean the GFF score and attribute data (`--clean` options)**

     The "clean" options perform the following operations:

     *   `--clean-score`: the data in the score column is cleaned up by replacing 'Anc_*'
         and blanks with '0's.

     The attribute field of the GFF can contain various semicolon-separated key-value
     pairs:

     *   `--clean-replace-attributes`: if one of these is a non-blank SGD then the Gene,
         Parent and Name values are updated to be the same as the SGD name.

     *   `--clean-exclude-attributes`: attributes called `kaks`, `kaks2` and `ncbi` are
         removed.

     If multiple features share the same SGD name then `--clean-replace-attributes` can
     result in them also sharing the same ID; to deal with this:

     *   `--clean-group-sgds`: update the ID attribute to group neighbouring lines that
         have the same SGD (see _"SGD grouping"_ below).

     A single `--clean` can be specified which performs all these operations
     automatically.

3.   **Detect duplicate SGDs (`--report-duplicates option`)**

     Report duplicate SGD names found in the input file.

     This option writes a list of the duplicates to a 'duplicates' file.

     It also reports the number of 'trivial' duplicates, i.e. lines having the same SGD
     because they are part of the same gene.

4.   **Resolve duplicate SGDs using a mapping file (`--resolve-duplicates option`)**

     Attempt to resolve duplicates by referring to a list of "best" genes given in a
     mapping file. For each duplicated name the resolution procedure is:

     *   Find mapping gene(s) with the same name

     *   For each mapping gene, keep duplicates which match chromosome, strand and which
         overlap with the start and end of the gene (see _"Overlap criteria"_ below). For
	 SGD groups the mapping gene must overlap the whole group for it to match; mapping
	 genes and duplicates which don't have matches are removed from the process.

     *   At the end of the matching procedure the duplication is resolved if there is one
         SGD (or SGD group) matched to one mapping gene. Otherwise the duplication remains
	 unresolved.

     When duplicates are resolved, the non-matching duplicates are discarded; otherwise by
     default all unresolved duplicates are kept. However if the `--discard-unresolved` option
     is also specified then all unresolved duplicates are removed before output; the
     `--insert-missing` option can then be used to add them back in.

     Note that the --discard-unresolved option cannot get rid of 'trivial' duplicates (i.e.
     lines having the same SGD because they are part of the same gene).

5.   **Add missing genes (`--insert-missing option`)**

     Adds genes from a list of "best" genes given in a mapping file which have names not
     found in the input GFF.

### Notes ###

#### SGD grouping ####

As part of setting the ID attribute of GFF lines, the "clean" option also attempts to group
neighbouring lines which have the same SGD name.

The ID attribute is updated to the form

        ID=CDS:<sgd_name>:<i>

where `<sgd_name>` is a gene or transcript name (e.g. `YEL0W`) and `<i>` is an integer
index which starts from 1. Groupings are indicated by subsequent lines having the same
`<sgd_name>` but monotonically increasing indices, for example:

        chr1   Test   CDS   34525   35262   0   -   0   ID=CDS:YEL0W:1;SGD=YEL0W
        chr1   Test   CDS   35823   37004   0   -   0   ID=CDS:YEL0W:2;SGD=YEL0W
        chr1   Test   CDS   38050   38120   0   -   0   ID=CDS:YEL0W:3;SGD=YEL0W
        chr1   Test   CDS   39195   39569   0   -   0   ID=CDS:YEL0W:4;SGD=YEL0W

When determining a grouping the program looks ahead from each line for subsequent lines
(up to five) which have the same SGD value. So groupings can also accommodate "breaks",
for example:

        chr1   Test   CDS   34525   35262   0   -   0   ID=CDS:YEL0W:1;SGD=YEL0W
        chr1   Test   CDS   35823   37004   0   -   0   ID=CDS:YEL0W:2;SGD=YEL0W
        chr1   Test   CDS   38050   38120   0   -   0   ID=CDS:YEL0X:1;SGD=YEL0X
        chr1   Test   CDS   39195   39569   0   -   0   ID=CDS:YEL0W:3;SGD=YEL0W

#### Mapping file format ####

The mapping file is a tab-delimited text file with lines of the form:

        name   chr   start   end    strand

`<name>` is used to match against the SGD names in the input GFF file.


#### Overlap criteria ####

Aside from matching chromosome and strand, one of the criteria for a mapping gene
to match a duplicate from the GFF file is that the two must overlap.

An overlap is counted as the duplicate from the GFF having start/end positions such
that it lies inside the start/end positions of the mapping gene extended by 1kb i.e.
between `start` - 1000 and `end` +  1000.

GFF3_Annotation_Extractor.py
----------------------------

GFF3_Annotation_Extractor.py takes gene feature data (for example the output from one
or more runs of the HTSeq-count program) and combines it with data about each feature's
parent gene from a GFF file.

By default the program takes a single tab-delimited input file where the first column
contains feature IDs, and appends data about the feature's parent gene.

In "htseq-count" mode, one or more `htseq-count` output files should be provided as
input, and the program will write out the data about the feature's parent gene appended
with the counts from each input file.

To generate the feature count files using `htseq-count` do e.g.:

        htseq-count --type=exon -i Parent <file>.gff <file>.sam

which returns counts of each exon against the name of that exon's parent.

GFF3_Annotation_Extractor will match up the exon parent with its parent gene and output
the counts against gene names.

### Usage ###

        GFF3_Annotation_Extractor.py OPTIONS <file>.gff FEATURE_DATA
	GFF3_Annotation_Extractor.py --htseq-count OPTIONS <file>.gff FEATURE_COUNTS [FEATURE_COUNTS2 ...]

### Options ###

    --version             show program's version number and exit
    -h, --help            show this help message and exit
    -o OUT_FILE           specify output file name
    -t FEATURE_TYPE, --type=FEATURE_TYPE
                          feature type contained in the process (default 'exon')
    --htseq-count         htseq-count mode: input is one or more output
                          FEATURE_COUNT files from the htseq-count program

### Output files ###

 * `<basename>_annot.txt`: the feature data annotated with data for each parent gene.

 * `<basename>_annot_stats.txt`: the counts of "ambiguous", "two_low_aQual" etc
    from each log (htseq-count mode only).

Set up and prerequisites
------------------------

The `GFFUtils` programs require the `TabFile.py` module from the FLS Bioinformatics
Core `genomics` repository.

Running tests
-------------

`GFFcleaner.py` has a set of unit tests built in; to run these do:

    GFFcleaner.py --test

To run the unit tests for the GFFFile.py module, do:

    python GFFFile.py