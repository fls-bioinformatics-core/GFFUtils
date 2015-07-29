#!/bin/sh
#
# Run test of GFF3_Annotation_Extractor.py with GTF input
# where we need to use -i gene_name
DATA_DIR=$(dirname $0)/data
if [ ! -d output ] ; then
  mkdir output
fi
GFF3_Annotation_Extractor.py \
 -o output/test_mm10_gencode_vM5.txt --htseq-count \
 -i gene_name \
 ${DATA_DIR}/mm10_gencode_vM5.gtf \
 ${DATA_DIR}/mm10_day6_s1_htseq_counts.txt
##
#
