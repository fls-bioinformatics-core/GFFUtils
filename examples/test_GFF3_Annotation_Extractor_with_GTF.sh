#!/bin/sh
#
# Run test of GFF3_Annotation_Extractor with GTF input
# Run test of GFF3_Annotation_Extractor
DATA_DIR=$(dirname $0)/data
if [ ! -d output ] ; then
  mkdir output
fi
GFF3_Annotation_Extractor \
 -o output/test_mm10.txt --htseq-count \
 ${DATA_DIR}/mm10.gtf \
 ${DATA_DIR}/mm10_htseq_counts.txt
##
#
