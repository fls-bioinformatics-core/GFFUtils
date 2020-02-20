#!/usr/bin/env python
#
#     gtf_extract.py: extract selected data from GTF file
#     Copyright (C) University of Manchester 2012 Peter Briggs
#

"""Utility program to extract selected data items from a GTF file.
"""

#######################################################################
# Module metadata
#######################################################################

from .. import get_version
__version__ = get_version()

#######################################################################
# Import modules
#######################################################################

import os
import sys
import logging
from argparse import ArgumentParser
from ..GFFFile import GFFIterator
from ..GFFFile import PRAGMA
from ..GFFFile import ANNOTATION
from ..GTFFile import GTFIterator

# Main program
#
def main():
    """Main program
    """
    # Command line parser
    p = ArgumentParser(description="Extract selected data items from a GTF file and "
                       "output in tab-delimited format. The program can also operate "
                       "on GFF files provided the --gff option is specified.")
    p.add_argument('gtf_file',metavar="GTF_FILE",
                   help="input GTF file to extract data items from")
    p.add_argument('-v','--version',action="version",version="%{prog}s "+__version__)
    p.add_argument('-f','--feature',action="store",dest="feature_type",default=None,
                   help="only extract data for lines where feature is FEATURE_TYPE")
    p.add_argument('--fields',action="store",dest="field_list",default=None,
                   help="comma-separated list of fields to output in tab-delimited "
                   "format for each line in the GTF, e.g. 'chrom,start,end'. Fields "
                   "can either be a GTF field name (i.e. 'chrom', 'source', "
                   "'feature', 'start', 'end', 'score', 'strand' and 'frame') or the "
                   "name of an attribute (e.g. 'gene_name', 'gene_id' etc). Data "
                   "items are output in the order they appear in FIELD_LIST. If a "
                   "field doesn't exist for a line then '.' will be output as the "
                   "value.")
    p.add_argument('-o',action="store",dest="outfile",default=None,
                   help="write output to OUTFILE (default is to write to stdout)")
    p.add_argument('--gff',action="store_true",dest="is_gff",default=False,
                   help="specify that the input file is GFF rather than GTF format")
    p.add_argument('-k','--keep-headers',action="store_true",dest="keep_header",
                   default=False,
                   help="copy headers from input file to output")
    args = p.parse_args()

    # Type of feature to extract data for
    feature_type = args.feature_type

    # Fields to report
    if args.field_list is None:
        field_list = None
    else:
        field_list = []
        for field in args.field_list.split(','):
            if field == "chr" or field == "chrom":
                # Allow 'chr' and 'chrom' to be aliases for 'seqname'
                field_list.append("seqname")
            else:
                field_list.append(field)

    # Input file type
    if args.is_gff:
        file_iterator =  GFFIterator
    else:
        file_iterator =  GTFIterator

    # Output stream
    if args.outfile is None:
        fp = sys.stdout
    else:
        fp = open(args.outfile,'w')

    # Null character (used when values are empty)
    null = '.'

    # Iterate through the file line-by-line
    for line in file_iterator(args.gtf_file):
        this_gene = None
        start = 0
        stop = 0
        if line.type == PRAGMA:
            if line[0].startswith('##gff-version') and not args.is_gff:
                logging.fatal("Input file is GFF not GTF? Rerun using --gff "
                              "option")
                sys.exit(1)
            if args.keep_header:
                fp.write("%s\n" % line)
        elif line.type == ANNOTATION:
            if feature_type is None or line['feature'] == feature_type:
                # Extract and report data
                if field_list is None:
                    fp.write("%s\n" % line)
                else:
                    out_line = []
                    for field in field_list:
                        try:
                            # Assume standard field
                            out_line.append(str(line[field]))
                        except KeyError:
                            # Not standard, try as an attribute name
                            try:
                                out_line.append(str(line['attributes'][field]))
                            except KeyError:
                                # Not an attribute either
                                out_line.append(str(null))
                    fp.write("%s\n" % '\t'.join(out_line))

    # Finished - close output file
    if args.outfile is not None:
        fp.close()

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s: %(message)s")
    main()


