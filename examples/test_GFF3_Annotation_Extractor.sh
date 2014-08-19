#!/bin/sh
#
# Run test of GFF3_Annotation_Extractor.py
if [ ! -d output ] ; then
  mkdir output
fi
GFF3_Annotation_Extractor.py -o output/test_dicty.out --htseq-count data/dicty.gff \
 data/dicty_htseq_counts.dimA.txt data/dicty_htseq_counts.AXA4.txt
##
#
