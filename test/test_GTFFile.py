#!/usr/bin/env python

import unittest
from io import StringIO
from GFFUtils.GTFFile import *

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

class TestGTFAttributes(unittest.TestCase):

    def setUp(self):
        #Example GTF data line
        self.gtf_line = """chr1	HAVANA	gene	11869	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"""

    def test_gtf_attributes(self):
        line = GTFDataLine(self.gtf_line)
        attributes = line['attributes']
        # Check the contents
        self.assertEqual(attributes['gene_id'],'ENSG00000223972.4')
        self.assertEqual(attributes['transcript_id'],'ENSG00000223972.4')
        self.assertEqual(attributes['gene_type'],'pseudogene')
        self.assertEqual(attributes['gene_status'],'KNOWN')
        self.assertEqual(attributes['gene_name'],'DDX11L1')
        self.assertEqual(attributes['transcript_type'],'pseudogene')
        self.assertEqual(attributes['transcript_status'],'KNOWN')
        self.assertEqual(attributes['transcript_name'],'DDX11L1')
        self.assertEqual(attributes['level'],'2')
        self.assertEqual(attributes['havana_gene'],'OTTHUMG00000000961.2')

    def test_gtf_contains(self):
        line = GTFDataLine(self.gtf_line)
        self.assertTrue('gene_id' in line['attributes'])
        self.assertTrue('transcript_id' in line['attributes'])
        self.assertTrue('gene_type' in line['attributes'])
        self.assertTrue('gene_status' in line['attributes'])
        self.assertTrue('gene_name' in line['attributes'])
        self.assertTrue('transcript_type' in line['attributes'])
        self.assertTrue('transcript_status' in line['attributes'])
        self.assertTrue('transcript_name' in line['attributes'])
        self.assertTrue('level' in line['attributes'])
        self.assertTrue('havana_gene' in line['attributes'])
        self.assertFalse('NONEXISTENT' in line['attributes'])

    def test_gtf_iteration(self):
        line = GTFDataLine(self.gtf_line)
        for attr in line['attributes']:
            self.assertNotEqual(line['attributes'][attr],None)

class TestGTFIterator(unittest.TestCase):
    """Basic tests for iterating through a GTF file
    """

    def setUp(self):
        # Example GTF file fragment
        self.fp = StringIO(
u"""##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)
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
        self.fp = StringIO(
u"""##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)
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
        # Check format and version
        self.assertEqual(gtf.format,'gtf')
        self.assertEqual(gtf.version,None)
        # There should be 6 lines of data
        self.assertEqual(len(gtf),6)
        # Check the feature names are correct
        #feature = ('chromosome','contig','gene','mRNA','exon','CDS')
        feature = ('gene','transcript','exon','exon','exon','transcript')
        for i in range(len(gtf)):
            self.assertEqual(feature[i],gtf[i]['feature'],
                             "Incorrect feature '%s' on data line %d" % (gtf[i]['feature'],i))
