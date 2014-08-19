#!/bin/sh
#
# Run test of GFF3_Annotation_Extractor.py with GTF input
# Run test of GFF3_Annotation_Extractor.py
if [ ! -d output ] ; then
  mkdir output
fi
GFF3_Annotation_Extractor.py -o output/test_mm10.txt --htseq-count data/mm10.gtf data/mm10_htseq_counts.txt
##
#
