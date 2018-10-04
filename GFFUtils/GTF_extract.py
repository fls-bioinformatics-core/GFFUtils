#!/bin/env python
#
#     GTF_extract.py: extract selected data from GTF file
#     Copyright (C) University of Manchester 2012 Peter Briggs
#
######################################################################
#
# GTF_extract.py
#
#######################################################################

"""GTF_extract

Utility program to extract selected data items from a GTF file.
"""

#######################################################################
# Module metadata
#######################################################################

from . import get_version
__version__ = get_version()

#######################################################################
# Import modules
#######################################################################

import os
import sys
import optparse
import GFFFile
import GTFFile

# Main program
#
def main():
    """Main program
    """
    # Command line parser
    p = optparse.OptionParser(usage="%prog OPTIONS <gft_file>",
                              version="%prog "+__version__,
                              description="Extract selected data items from a GTF file and "
                              "output in tab-delimited format. The program can also operate "
                              "on GFF files provided the --gff option is specified.")
    p.add_option('-f','--feature',action="store",dest="feature_type",default=None,
                 help="only extract data for lines where feature is FEATURE_TYPE")
    p.add_option('--fields',action="store",dest="field_list",default=None,
                 help="comma-separated list of fields to output in tab-delimited format "
                 "for each line in the GTF, e.g. 'chrom,start,end'. Fields can either be a "
                 "GTF field name (i.e. 'chrom', 'source', 'feature', 'start', 'end', 'score', "
                 "'strand' and 'frame') or the name of an attribute (e.g. 'gene_name', "
                 "'gene_id' etc). Data items are output in the order they appear in "
                 "FIELD_LIST. If a field doesn't exist for a line then '.' will be output "
                 "as the value.")
    p.add_option('-o',action="store",dest="outfile",default=None,
                 help="write output to OUTFILE (default is to write to stdout)")
    p.add_option('--gff',action="store_true",dest="is_gff",default=False,
                 help="specify that the input file is GFF rather than GTF format")
    p.add_option('-k','--keep-headers',action="store_true",dest="keep_header",default=False,
                 help="copy headers from input file to output")

    opts,args = p.parse_args()

    # Check number of arguments
    if len(args) != 1:
        p.error("Expected single argument (GTF file)")

    # Type of feature to extract data for
    feature_type = opts.feature_type

    # Fields to report
    if opts.field_list is None:
        field_list = None
    else:
        field_list = []
        for field in opts.field_list.split(','):
            if field == "chr" or field == "chrom":
                # Allow 'chr' and 'chrom' to be aliases for 'seqname'
                field_list.append("seqname")
            else:
                field_list.append(field)

    # Input file type
    if opts.is_gff:
        file_iterator =  GFFFile.GFFIterator
    else:
        file_iterator =  GTFFile.GTFIterator

    # Output stream
    if opts.outfile is None:
        fp = sys.stdout
    else:
        fp = open(opts.outfile,'w')

    # Null character (used when values are empty)
    null = '.'

    # Iterate through the file line-by-line
    for line in file_iterator(args[0]):
        this_gene = None
        start = 0
        stop = 0
        if line.type == GFFFile.PRAGMA:
            if line[0].startswith('##gff-version') and not opts.is_gff:
                sys.stderr.write("Input file is GFF not GTF? Rerun using --gff option\n")
                sys.exit(1)
            if opts.keep_header:
                fp.write("%s\n" % line)
        elif line.type == GFFFile.ANNOTATION:
            if feature_type is None or line['feature'] == feature_type:
                # Extract and report data
                if field_list is None:
                    print "%s" % line
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
    if opts.outfile is not None:
        fp.close()

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s: %(message)s")
    main()


