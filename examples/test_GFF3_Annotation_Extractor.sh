#!/bin/sh
#
# Run test of GFF3_Annotation_Extractor
DATA_DIR=$(dirname $0)/data
if [ ! -d output ] ; then
  mkdir output
fi
GFF3_Annotation_Extractor \
 -o output/test_dicty.out --htseq-count \
 ${DATA_DIR}/dicty.gff \
 ${DATA_DIR}/dicty_htseq_counts.dimA.txt \
 ${DATA_DIR}/dicty_htseq_counts.AXA4.txt
##
#
