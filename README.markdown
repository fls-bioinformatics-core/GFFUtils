GFFUtils
========

GFFUtils package provides the following utilities for working with GFF files:

 * `GFFcleaner.py`: performs various "cleaning" manipulations on a GFF file
 * `GFF_HTSeq_Annotator.py`: combine and annotate HTSeq-count output with data from
   the source GFF file

GFFcleaner.py
-------------

The GFFcleaner can perform various manipulations on a GFF file to "clean" it.

### Usage ###

        GFFcleaner.py [OPTIONS] <file>.gff

### Options ###

 *  `-o OUTPUT_GFF`

    Explicitly set the output file name (otherwise defaults to `<file>_clean.gff`)

 *  `--prepend=<str>`

    String to prepend to seqname in first column

 *  `--clean`

    Perform the 'cleaning' manipulations on data

 *  `--report-duplicates`

    Report duplicate SGD names

 *  `--resolve-duplicates=<mapping-file>`

    Resolve duplicates by matching against 'best' genes
    in `<mapping-file>`. Non-matching genes are discarded.

 *  `--discard-unresolved`

    Also discard any unresolved duplicates

 *  `--insert-missing[=<mapping-file>]`

    Insert genes from <mapping-file> with SGD names
    that don't appear in the input GFF

    (If `<mapping-file>` is specified with the
    `--resolve-duplicates` option then that will be
    used by default.)

 *  `--add-exon-ids`

    For exon features without an ID attribute, construct
    and insert an ID of the form `exon:<Parent>:<n>`
    (where n is a unique number).

 *  `--add-missing-ids`

     For features without an ID attribute, construct and
     insert a generated ID of the form
     `<feature>:<Parent>:<n>` (where n is a unique number).

 *  `--no-percent-encoding`

     Convert encoded attributes to the correct characters
     in the output GFF. **WARNING this may result in a non-
     cannonical GFF that can't be read correctly by this or
     other programs.**


 *  `--debug` Print debugging information

 *  `--test`  Run unit tests

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

2.   **Clean the GFF score and attribute data (`--clean option`)**

     The clean option does a number of things:

     *   The data in the score column is cleaned up by replacing 'Anc_*' and blanks with '0's.

     The attribute field of the GFF can contain various semicolon-separated key-value
     pairs:

     *   If one of these is a non-blank SGD then the Gene, Parent and Name values are
         updated to be the same as the SGD name; ID is also updated to group neighbouring
	 lines with the same SGD (see _"SGD grouping"_ below).

     *   Attributes called `kaks`, `kaks2` and `ncbi` are removed.

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

GFF_HTSeq_Annotator.py
----------------------

GFF_HTSeq_Annotator.py takes the output from one or more runs of the HTSeq-count program
and combines it with data from the source GFF file.

HTSeq-count should have been run using e.g.:

        htseq-count --type=exon -i Parent <file>.gff <file>.sam

which returns counts of each exon against the name of that exon's parent.

GFF_HTSeq_Annotator will match up the exon parent with its parent gene and output the
counts against gene names.

### Usage ###

        GFF_HTSeq_Annotator.py <file>.gff <htseq-log1> [<htseq-log2> ...]

Alternatively you can use wildcards to specify multiple HTSeq log files e.g. `*_HTSeq.txt`.

### Options ###

    --version             show program's version number and exit
    -h, --help            show this help message and exit
    -o OUT_FILE           specify output file name
    -t FEATURE_TYPE, --type=FEATURE_TYPE
                          feature type to process (default 'exon')

### Output files ###

 * `<basename>_htseq_counts.txt`: the counts from each log for each gene name.

 * `<basename>_htseq_counts_stats.txt`: the counts of "ambiguous", "two_low_aQual" etc
    from each log.

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