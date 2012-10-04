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

import version
__version__ = version.__version__

#######################################################################
# Import modules
#######################################################################

import os
import sys
import optparse
import GFFFile

#######################################################################
# Classes
#######################################################################

class GTFDataLine(GFFFile.GFFDataLine):
    """Data line specific to GTF files

    Subclass of GFFDataLine which adds two properties:

    gene_id: returns the gene_id value stored in the attributes
    transcript_id: returns the transcript_id value stored in the attributes
    """

    def __init__(self,line=None,column_names=GFFFile.GFF_COLUMNS,lineno=None,delimiter='\t',
                 gff_line_type=None):
        GFFFile.GFFDataLine.__init__(self,line=line,column_names=column_names,
                                     lineno=lineno,delimiter=delimiter,gff_line_type=gff_line_type)
        self.__attributes = {}

    def attribute(self,name):
        """
        """
        if name in self.__attributes:
            return self.__attributes[name]
        else:
            for attr in self['attributes'].nokeys():
                if attr.startswith("%s " % name):
                    self.__attributes[name] = attr.split(' ')[1].strip('"')
                    return self.__attributes[name]
        return None

#######################################################################
# Tests
#######################################################################

import unittest
import cStringIO

class TestGTFDataLine(unittest.TestCase):

    def setUp(self):
        # Example GTF data line
        self.gtf_line = """chr1	HAVANA	gene	11869	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"""

    def test_gtf_data_line(self):
        line = GTFDataLine(self.gtf_line)
        # Check basic data items
        self.assertEqual("chr1",line['seqname'])
        self.assertEqual("HAVANA",line['source'])
        self.assertEqual("gene",line['feature'])
        self.assertEqual(11869,line['start'])
        self.assertEqual(14412,line['end'])
        self.assertEqual(".",line['score'])
        self.assertEqual("+",line['strand'])
        self.assertEqual(".",line['frame'])
        # Check attributes
        self.assertEqual("ENSG00000223972.4",line.attribute('gene_id'))
        self.assertEqual("ENSG00000223972.4",line.attribute('transcript_id'))
        self.assertEqual("pseudogene",line.attribute('gene_type'))
        self.assertEqual("KNOWN",line.attribute('gene_status'))
        self.assertEqual("DDX11L1",line.attribute('gene_name'))
        self.assertEqual("pseudogene",line.attribute('transcript_type'))
        self.assertEqual("KNOWN",line.attribute('transcript_status'))
        self.assertEqual("DDX11L1",line.attribute('transcript_name'))
        self.assertEqual("2",line.attribute('level'))
        self.assertEqual("OTTHUMG00000000961.2",line.attribute('havana_gene'))

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Command line parser
    p = optparse.OptionParser(usage="%prog OPTIONS <gft_file>",
                              version="%prog "+__version__,
                              description="Extract selected data items from a GTF file and "
                              "output in tab-delimited format.")
    p.add_option('-f','--feature',action="store",dest="feature_type",default=None,
                 help="only extract data for lines where feature is FEATURE_TYPE")
    p.add_option('--fields',action="store",dest="field_list",default=None,
                 help="comma-separated list of fields to output in tab-delimited format "
                 "for each line in the GTF, e.g. 'chrom,start,end'. Fields can either be a "
                 "GTF field name (i.e. 'chrom', 'source', 'feature', 'start', 'end', 'score', "
                 "'strand' and 'frame') or the name of an attribute (e.g. 'gene_name', "
                 "'gene_id' etc). Data items are output in the order they appear in "
                 "FIELD_LIST.")
    p.add_option('-o',action="store",dest="outfile",default=None,
                 help="write output to OUTFILE (default is to write to stdout)")
    p.add_option('--test',action="store_true",dest="run_tests",default=False,
                 help="run unit tests (developers only)")

    opts,args = p.parse_args()

    # Check for unit testing
    if opts.run_tests:
        print "Running unit tests"
        suite = unittest.TestSuite(unittest.TestLoader().\
                                       discover(os.path.dirname(sys.argv[0]), \
                                                    pattern=os.path.basename(sys.argv[0])))
        unittest.TextTestRunner(verbosity=2).run(suite)
        print "Tests finished"
        sys.exit()

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

    # Output stream
    if opts.outfile is None:
        fp = sys.stdout
    else:
        fp = open(opts.outfile,'w')

    # Iterate through the GTF file line-by-line
    for line in GFFFile.GFFIterator(args[0],gffdataline=GTFDataLine):
        this_gene = None
        start = 0
        stop = 0
        if line.type == GFFFile.ANNOTATION:
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
                            out_line.append(str(line.attribute(field)))
                    fp.write("%s\n" % '\t'.join(out_line))

    # Finished - close output file
    if opts.outfile is not None:
        fp.close()
                    
        
