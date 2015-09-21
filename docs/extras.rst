Extra scripts and utilities
===========================

run_cleanup.sh
--------------

``run_cleanup.sh`` is a shell script which automatically runs
``GFFcleaner`` using all the steps outlined in the _Usage recipe_.

Usage::

    run_cleanup.sh <gff_file> [<mapping_file>]

Note that duplicate resolution and missing gene insertion cannot
be performed unless a mapping file is supplied.
