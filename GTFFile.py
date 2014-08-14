#!/bin/env python
#
#     GFTFile.py: classes for reading and manipulating data in GTF files
#     Copyright (C) University of Manchester 2012-4 Peter Briggs
#
########################################################################
#
# GFTFile.py
#
#########################################################################

import version
__version__ = version.__version__

"""GTFFile

Classes for reading data from GTF (Gene Transfer Format/General Transfer
Format) files.

According to the Ensembl documentation, GTF is the same as GFF version 2:

http://www.ensembl.org/info/website/upload/gff.html

however specific documentation for GTF2 can also be found at

http://mblab.wustl.edu/GTF2.html

Classes
-------

There are a number of GTF-specific classes:

 * GTFIterator: line-by-line iteration through a GTF
 * GTFFile: read data from GTF into memory so it can be easily interrogated
 * GTFDataLine: get data from a single line from a GTF file

These classes are built on top of the GFF handling classes.

Usage examples
--------------

Start by reading the data from a GTF file into a GTFFile object:

>>> gtf = GTFFile('my.gtf')

To iterate over all lines in the GTFFile object and print the
feature type:

>>> for line in gtf:
>>>    print line['feature']

(To iterate over all lines in the GTF without caching in memory,
printing feature types for annotation records:

>>> for line in GTFIterator('my.gtf'):
>>>    if line.type == ANNOTATION:
>>>       print line['feature']

)

The 'attributes' data of each annotation record are automatically converted
to a GFFAttributes object, which allows the values of named attributes to be
referenced directly. For example, if the attributes string is:

"ID=1690892742066571889;SGD=YEL026W;Gene=Sbay_5.43;Parent=Sbay_5.43"

then the value of the 'Parent' attribute can be obtained using

>>> print line['attributes']['Parent']

To iterate over all lines and extract 'ID' field from attributes:

>>> for line in gff:
>>>    if 'ID' in line['attributes']['ID']:
>>>       print line['attributes']['ID']

To iterate over all lines and print just the 'name' part of the
'ID' attribute:

>>> for line in gff:
>>>    if 'ID' in line['attributes']
>>>       print GFFID(line['attributes']['ID']).name

"""

import GFFFile

#######################################################################
# Classes
#######################################################################

class GTFDataLine(GFFFile.GFFDataLine):
    """Data line specific to GTF files

    Subclass of GFFDataLine specifically for handling GTF data.

    """

    def __init__(self,line=None,column_names=GFFFile.GFF_COLUMNS,lineno=None,delimiter='\t',
                 gff_line_type=None):
        GFFFile.GFFDataLine.__init__(self,line=line,column_names=column_names,
                                     lineno=lineno,delimiter=delimiter,gff_line_type=gff_line_type)
        self['attributes'] = GTFAttributes(str(self['attributes']))

class GTFAttributes:
    def __init__(self,attribute_data=None):
        self.__attributes = GFFFile.OrderedDictionary()
        for attr in GFFFile.GFFAttributes(attribute_data=attribute_data).nokeys():
            key = attr.split(' ')[0]
            self.__attributes[key] = ' '.join(attr.split(' ')[1:]).strip('"')
    def __getitem__(self,name):
        try:
            return self.__attributes[name]
        except KeyError:
            return None
    def __contains__(self,name):
        return self[name] is not None

class GTFFile(GFFFile.GFFFile):
    """Class for handling GTF files in-memory

    Subclass of GFFFile which uses a GTFDataLine to store the
    annotation lines.

    Data from the file can then be extracted and modified using the
    methods of the GFFFile superclass (and its TabFile superclass)
    and the GTFDataLine.

    GTF is alledgedly the same as GFF version 2. See
    http://www.sanger.ac.uk/resources/software/gff/spec.html
    for the GFF specification.

    """
    def __init__(self,gtf_file,fp=None,**args):
        args['gffdataline'] = GTFDataLine
        GFFFile.GFFFile.__init__(self,gtf_file,fp=fp,**args)

class GTFIterator(GFFFile.GFFIterator):
    def __init__(self,gtf_file=None,fp=None,**args):
        args['gffdataline'] = GTFDataLine
        GFFFile.GFFIterator.__init__(self,gff_file=gtf_file,fp=fp,**args)

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
        self.assertEqual("ENSG00000223972.4",line['attributes']['gene_id'])
        self.assertEqual("ENSG00000223972.4",line['attributes']['transcript_id'])
        self.assertEqual("pseudogene",line['attributes']['gene_type'])
        self.assertEqual("KNOWN",line['attributes']['gene_status'])
        self.assertEqual("DDX11L1",line['attributes']['gene_name'])
        self.assertEqual("pseudogene",line['attributes']['transcript_type'])
        self.assertEqual("KNOWN",line['attributes']['transcript_status'])
        self.assertEqual("DDX11L1",line['attributes']['transcript_name'])
        self.assertEqual("2",line['attributes']['level'])
        self.assertEqual("OTTHUMG00000000961.2",line['attributes']['havana_gene'])
        self.assertEqual(None,line['attributes']['missing'])

class TestGTFIterator(unittest.TestCase):
    """Basic tests for iterating through a GTF file
    """

    def setUp(self):
        # Example GTF file fragment
        self.fp = cStringIO.StringIO(
"""##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)
##provider: GENCODE
##contact: gencode@sanger.ac.uk
##format: gtf
##date: 2013-12-05
chr1	HAVANA	gene	11869	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 1;  exon_id "ENSE00002234944.1";  level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 2;  exon_id "ENSE00003582793.1";  level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 3;  exon_id "ENSE00002312635.1";  level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	ENSEMBL	transcript	11872	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000515242.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1-201"; level 3; havana_gene "OTTHUMG00000000961.2";
""")

    def test_gtf_iterator(self):
        """Test iteration over a file-like object
        """
        # Count number of lines, and number of each type
        nlines = 0
        npragma = 0
        ncomment = 0
        nannotation = 0
        # Iterate through GTF file fragment
        for line in GTFIterator(fp=self.fp):
            nlines += 1
            self.assertNotEqual(line.type,None)
            self.assertEqual(line.lineno(),nlines)
            if line.type == GFFFile.PRAGMA: npragma += 1
            if line.type == GFFFile.COMMENT: ncomment += 1
            if line.type == GFFFile.ANNOTATION: nannotation += 1
        # Check counts
        self.assertEqual(nlines,11)
        self.assertEqual(npragma,5)
        self.assertEqual(ncomment,0)
        self.assertEqual(nannotation,6)

class TestGTFFile(unittest.TestCase):
    """Basic unit tests for the GTFFile class
    """

    def setUp(self):
        # Example GTF file fragment
        self.fp = cStringIO.StringIO(
"""##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)
##provider: GENCODE
##contact: gencode@sanger.ac.uk
##format: gtf
##date: 2013-12-05
chr1	HAVANA	gene	11869	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 1;  exon_id "ENSE00002234944.1";  level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 2;  exon_id "ENSE00003582793.1";  level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 3;  exon_id "ENSE00002312635.1";  level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1	ENSEMBL	transcript	11872	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000515242.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1-201"; level 3; havana_gene "OTTHUMG00000000961.2";
""")

    def test_read_in_gtf(self):
        """Test that the GTF data can be read in
        """
        gtf = GTFFile("test.gtf",self.fp)
        # There should be 6 lines of data
        self.assertEqual(len(gtf),6)
        # Check the feature names are correct
        #feature = ('chromosome','contig','gene','mRNA','exon','CDS')
        feature = ('gene','transcript','exon','exon','exon','transcript')
        for i in range(len(gtf)):
            self.assertEqual(feature[i],gtf[i]['feature'],
                             "Incorrect feature '%s' on data line %d" % (gtf[i]['feature'],i))
