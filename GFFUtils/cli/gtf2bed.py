#!/usr/bin/env python
#
#     gtf2bed: convert GTF contents to BED format
#     Copyright (C) University of Manchester 2012,2020 Peter Briggs
#
import sys
import logging
from argparse import ArgumentParser
from ..GTFFile import GTFIterator
from ..GFFFile import ANNOTATION
from .. import get_version

def main():
    """
    gtf2bed: convert GTF file to BED format

    """
    p = ArgumentParser(description="Convert GTF file to BED format")
    p.add_argument('gtf_in',metavar="FILE.gtf",
                   help="GTF file to convert")
    p.add_argument('-v','--version',action='version',
                   version=get_version())
    p.add_argument('-o',action="store",dest="outfile",default=None,
                   help="write output to OUTFILE (default is to write "
                   "to stdout)")
    args = p.parse_args()
    # Output stream
    if args.outfile is None:
        fp = sys.stdout
    else:
        fp = open(args.outfile,'wt')
    # Iterate over GTF
    for line in GTFIterator(args.gtf_in):
        this_gene = None
        start = 0
        stop = 0
        attributes = line['attributes']
        if line.type == ANNOTATION:
            if line['feature'] == 'gene':
                # Encountered gene feature
                if this_gene != attributes['gene_name']:
                    this_gene = attributes['gene_name']
                    start = line['start']
                    stop = line['end']
                    fp.write("%s\t%s\t%d\t%d\t%s\t%s\n" %
                             (attributes['gene_name'],
                              line['seqname'],
                              line['start'],
                              line['end'],
                              line['strand'],
                              line['feature']))
            else:
                # Non-gene feature
                if this_gene == attributes['gene_name']:
                    print_details = False
                    if line['start'] < start:
                        logging.warn("WARNING Start is before gene start")
                        print_details = True
                    elif line['end'] > stop:
                        logging.warn("WARNING End is after gene end "
                                     "(%s > %s)" % (line['end'],stop))
                        print_details = True
                    if line['end'] < start:
                        logging.warn("WARNING End is before gene start")
                        print_details = True
                    elif line['start'] > stop:
                        logging.warn("WARNING Start is after gene end")
                        print_details = True
                    if print_details:
                        logging.warn("%s\t%s\t%s\t%s\t%s\t%s" %
                                     (attributes['gene_name'],
                                      line['seqname'],
                                      line['start'],
                                      line['end'],
                                      line['strand'],
                                      line['feature']))
    fp.close()
        
if __name__ == "__main__":
    main()
