``gtf2bed``: convert GTF contents to BED format
===============================================

Overview
--------

``gtf2bed`` converts the contents of a GTF file to BED format, printing
a single line for each ``gene`` entry in the input GTF.

Usage
-----

General usage syntax::

    gtf2bed FILE.gtf

Output
------

The program prints the BED file contents directly to stdout, for example::

    Gnai3	3	108107280	108146146	-	gene
    Pbsn	X	77837901	77853623	-	gene
    Cdc45	16	18780447	18811987	-	gene
    H19	7	142575529	142578143	-	gene
    Scml2	X	161117193	161258213	+	gene

To sent this to a file use::

    gtf2bed FILE.gtf > FILE.bed

If any problems are encountered then the program will print warning
messages to stderr.
