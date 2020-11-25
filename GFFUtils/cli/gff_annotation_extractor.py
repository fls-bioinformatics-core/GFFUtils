#!/usr/bin/env python
#
#     gff3_annotation_extractor.py: annotate feature counts with data from GFF
#     Copyright (C) University of Manchester 2012,2020 Peter Briggs, Leo Zeef
#

"""
Annotate gene feature data (for example the output from one or more runs of the
HTSeq-count program) by combining it with data about each feature's parent gene,
taken from a GFF file.

By default the program takes a single tab-delimited input file where the first column
contains feature IDs, and appends data about the feature's parent gene.

In "htseq-count" mode, one or more `htseq-count` output files should be provided as
input, and the program will write out the data about the feature's parent gene appended
with the counts from each input file.

To generate the feature count files using `htseq-count` do e.g.:

        htseq-count --type=exon -i Parent <file>.gff <file>.sam

which returns counts of each exon against the name of that exon's parent.

Multiple parents
----------------

It's possible for features in a GFF file to have multiple parents.

In this case the output from htseq-count will reproduce the 'Parent'
attribute verbatim, e.g. AF2312,AB2812,abc-3. However
GFF3_Annotation_Extractor will be unable to determine the parent genes
in this case, and so will issue a warning and continue.
"""

#######################################################################
# Module metadata
#######################################################################

from .. import get_version
__version__ = get_version()

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import os
import glob
import logging
from argparse import ArgumentParser
from ..GFFFile import GFFFile
from ..GTFFile import GTFFile
from ..annotation import GFFAnnotationLookup
from ..annotation import annotate_htseq_count_data
from ..annotation import annotate_feature_data

# Main program
#
def main():
    """Main program
    """
    # Process command line
    p = ArgumentParser(description="Annotate feature count data with information "
                       "from a GFF or GTF file. Default mode is to take a single "
                       "tab-delimited FEATURE_DATA input file where the first column "
                       "consists of feature IDs from the input GFF_FILE; in this "
                       "mode each line of FEATURE_DATA will be appended with data "
                       "about the 'parent feature' and 'parent gene' matching the "
                       "feature ID. In --htseq-count mode input consists of one or "
                       "more FEATURE_COUNTS files generated using htseq-count (e.g. "
                       "'htseq-count -q -t exon -i Parent gff_file sam_file'). The "
                       "annotator looks up the parent genes of each feature and "
                       "outputs this information against the feature counts "
                       "(in <GFF_FILE>_annot.txt) plus the totals assigned, not "
                       "counted etc (in <GFF_FILE>_annot_stats.txt).")
    p.add_argument('-v','--version',action='version',version=__version__)
    p.add_argument('gff_file',metavar="GFF_FILE",
                   help="GFF or GTF file to get annotation data from")
    p.add_argument('feature_files',metavar="FEATURE_FILE",nargs="+",
                   help="feature data to annotate; should be either "
                   "tab-delimited data, or output from htseq-count (if "
                   "--htseq-count is also specified)")
    p.add_argument('-o',action="store",dest="out_file",default=None,
                   help="specify output file name")
    p.add_argument('-t','--type',action="store",dest="feature_type",
                   default=None,
                   help="restrict feature records to this type when "
                   "matching features from input count files; if used in "
                   "conjunction with --htseq-count then should be the "
                   "same as that specified when running htseq-count "
                   "(default: include all feature records)")
    p.add_argument('-i','--id-attribute',action="store",dest="id_attribute",
                   default=None,
                   help="explicitly specify the name of the attribute to get the "
                   "feature IDs from (defaults to 'ID' for GFF input, 'gene_id' for "
                   "GTF input)")
    p.add_argument('--htseq-count',action="store_true",dest="htseq_count",
                   default=False,
                   help="htseq-count mode: input is one or more FEATURE_FILEs "
                   "output from htseq-count")
    args = p.parse_args()

    # Determine what mode to operate in
    htseq_count_mode = args.htseq_count

    # Input GFF file
    gff_file = args.gff_file
    if not os.path.exists(gff_file):
        p.error("Input GFF/GTF file %s not found" % gff_file)

    # Check for wildcards in feature data file names, to emulate linux shell globbing
    # on platforms such as Windows which don't have this built in
    feature_data_files = []
    for arg in args.feature_files:
        for filen in glob.iglob(arg):
            if not os.path.exists(filen):
                p.error("File '%s' not found" % filen)
            feature_data_files.append(filen)
    if not feature_data_files:
        p.error("No input feature data files found")

    # Final check on number of input files
    if not htseq_count_mode and len(feature_data_files) > 1:  
        p.error("Expected GFF/GTF file and a single feature data file")

    # Feature type being considered
    feature_type = args.feature_type

    # Output file
    if args.out_file:
        out_file = args.out_file
    else:
        out_file = os.path.splitext(os.path.basename(gff_file))[0] + "_annot.txt"

    # Process GFF/GTF data
    print("Reading data from %s" % gff_file)
    if gff_file.endswith('.gtf'):
        gff = GTFFile(gff_file)
    else:
        gff = GFFFile(gff_file)
    feature_format = gff.format.upper()

    # Build lookup
    print("Creating lookup for %s" % feature_format)
    feature_lookup = GFFAnnotationLookup(gff,
                                         id_attr=args.id_attribute,
                                         feature_type=feature_type)

    # Annotate input data
    if htseq_count_mode:
        # HTSeq-count mode
        annotate_htseq_count_data(feature_lookup,
                                  feature_data_files,
                                  out_file)
    else:
        # Standard mode
        annotate_feature_data(feature_lookup,
                              feature_data_files[0],
                              out_file)

def GFF3_Annotation_Extractor():
    """
    Deprecated frontend
    """
    logging.warning("'GFF3_Annotation_Extractor' is deprecated and will "
                    "be removed in a future release")
    logging.warning("Please use 'gff_annotation_extractor'")
    main()

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    main()
