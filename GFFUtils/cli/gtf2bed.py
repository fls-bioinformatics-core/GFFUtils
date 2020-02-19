#!/usr/bin/env python
#
#     gtf2bed: convert GTF contents to BED format
#     Copyright (C) University of Manchester 2012 Peter Briggs
#
import sys
import optparse
from ..GTFFile import GTFIterator
from ..GFFFile import ANNOTATION

def warning(msg):
    """
    Print message to stderr
    """
    sys.stderr.write("%s\n" % msg)

def main():
    """
    gtf2bed: convert GTF file to BED format

    """
    p = optparse.OptionParser(usage="%prog FILE.gtf",
                              description="Convert GTF file to BED format")
    p.add_option('-o',action="store",dest="outfile",default=None,
                 help="write output to OUTFILE (default is to write to stdout)")
    opts,args = p.parse_args()
    if len(args) != 1:
        p.error("Expected single argument (GTF file)")
    # Output stream
    if opts.outfile is None:
        fp = sys.stdout
    else:
        fp = open(opts.outfile,'w')
    # Iterate over GTF
    for line in GTFIterator(args[0]):
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
                        warning("WARNING Start is before gene start")
                        print_details = True
                    elif line['end'] > stop:
                        warning("WARNING End is after gene end (%s > %s)" %
                                (line['end'],stop))
                        print_details = True
                    if line['end'] < start:
                        warning("WARNING End is before gene start")
                        print_details = True
                    elif line['start'] > stop:
                        warning("WARNING Start is after gene end")
                        print_details = True
                    if print_details:
                        warning("%s\t%s\t%s\t%s\t%s\t%s" %
                                (attributes['gene_name'],
                                 line['seqname'],
                                 line['start'],
                                 line['end'],
                                 line['strand'],
                                 line['feature']))
    fp.close()
        
if __name__ == "__main__":
    main()
