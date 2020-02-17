#!/usr/bin/env python

import unittest
from io import StringIO
from GFFUtils.GFFFile import *

class TestGFFIterator(unittest.TestCase):
    """Basic tests for iterating through a GFF file
    """

    def setUp(self):
        # Example GFF file fragment
        self.fp = StringIO(
"""##gff-version 3
# generated: Wed Feb 21 12:01:58 2012
DDB0123458	Sequencing Center	chromosome	1	4923596	.	+	.	ID=DDB0232428;Name=1
DDB0232428	Sequencing Center	contig	101	174493	.	+	.	ID=DDB0232440;Parent=DDB0232428;Name=DDB0232440;description=Contig generated from contig finding genome version 2.5;Dbxref=Contig GI Number:90970918,Accession Number:AAFI02000001,SeqID for Genbank:DDB0232440.02
DDB0232428	.	gene	1890	3287	.	+	.	ID=DDB_G0267178;Name=DDB_G0267178_RTE;description=ORF2 protein fragment of DIRS1 retrotransposon%3B refer to Genbank M11339 for full-length element
DDB0232428	Sequencing Center	mRNA	1890	3287	.	+	.	ID=DDB0216437;Parent=DDB_G0267178;Name=DDB0216437;description=JC1V2_0_00003: Obtained from the Dictyostelium Genome Consortium at The Wellcome Trust Sanger Institute;translation_start=1;Dbxref=Protein Accession Version:EAL73826.1,Inparanoid V. 5.1:DDB0216437,Protein Accession Number:EAL73826.1,Protein GI Number:60475899,UniProt:Q55H43,Genome V. 2.0 ID:JC1V2_0_00003
DDB0232428	Sequencing Center	exon	1890	3287	.	+	.	Parent=DDB0216437
DDB0232428	Sequencing Center	CDS	1890	3287	.	+	.	Parent=DDB0216437
""")

    def test_gff_iterator(self):
        """Test iteration over a file-like object
        """
        # Count number of lines, and number of each type
        nlines = 0
        npragma = 0
        ncomment = 0
        nannotation = 0
        # Iterate through GFF file fragment
        for line in GFFIterator(fp=self.fp):
            nlines += 1
            self.assertNotEqual(line.type,None)
            self.assertEqual(line.lineno(),nlines)
            if line.type == PRAGMA: npragma += 1
            if line.type == COMMENT: ncomment += 1
            if line.type == ANNOTATION: nannotation += 1
        # Check counts
        self.assertEqual(nlines,8)
        self.assertEqual(npragma,1)
        self.assertEqual(ncomment,1)
        self.assertEqual(nannotation,6)

class TestGFFFile(unittest.TestCase):
    """Basic unit tests for the GFFFile class
    """

    def setUp(self):
        # Example GFF file fragment
        self.fp = StringIO(
"""##gff-version 3
# generated: Wed Feb 21 12:01:58 2012
DDB0123458	Sequencing Center	chromosome	1	4923596	.	+	.	ID=DDB0232428;Name=1
DDB0232428	Sequencing Center	contig	101	174493	.	+	.	ID=DDB0232440;Parent=DDB0232428;Name=DDB0232440;description=Contig generated from contig finding genome version 2.5;Dbxref=Contig GI Number:90970918,Accession Number:AAFI02000001,SeqID for Genbank:DDB0232440.02
DDB0232428	.	gene	1890	3287	.	+	.	ID=DDB_G0267178;Name=DDB_G0267178_RTE;description=ORF2 protein fragment of DIRS1 retrotransposon%3B refer to Genbank M11339 for full-length element
DDB0232428	Sequencing Center	mRNA	1890	3287	.	+	.	ID=DDB0216437;Parent=DDB_G0267178;Name=DDB0216437;description=JC1V2_0_00003: Obtained from the Dictyostelium Genome Consortium at The Wellcome Trust Sanger Institute;translation_start=1;Dbxref=Protein Accession Version:EAL73826.1,Inparanoid V. 5.1:DDB0216437,Protein Accession Number:EAL73826.1,Protein GI Number:60475899,UniProt:Q55H43,Genome V. 2.0 ID:JC1V2_0_00003
DDB0232428	Sequencing Center	exon	1890	3287	.	+	.	Parent=DDB0216437
DDB0232428	Sequencing Center	CDS	1890	3287	.	+	.	Parent=DDB0216437
""")

    def test_read_in_gff(self):
        """Test that the GFF data can be read in
        """
        gff = GFFFile("test.gff",self.fp)
        # Check format and version
        self.assertEqual(gff.format,'gff')
        self.assertEqual(gff.version,'3')
        # There should be 6 lines of data
        self.assertEqual(len(gff),6)
        # Check the feature names are correct
        feature = ('chromosome','contig','gene','mRNA','exon','CDS')
        for i in range(len(gff)):
            self.assertEqual(feature[i],gff[i]['feature'],
                             "Incorrect feature '%s' on data line %d" % (gff[i]['feature'],i))

class TestGFFAttributes(unittest.TestCase):
    """Unit tests for GFFAttributes class
    """

    def test_extract_attributes(self):
        """Test that attribute data can be extracted
        """
        attributes = "ID=DDB0232440;Parent=DDB0232428;Name=DDB0232440;description=Contig generated from contig finding genome version 2.5;Dbxref=Contig GI Number:90970918,Accession Number:AAFI02000001,SeqID for Genbank:DDB0232440.02"
        attr = GFFAttributes(attributes)
        self.assertTrue('ID' in attr)
        self.assertEqual(attr['ID'],'DDB0232440')
        self.assertTrue('Parent' in attr)
        self.assertEqual(attr['Parent'],'DDB0232428')
        self.assertTrue('Name' in attr)
        self.assertEqual(attr['Name'],'DDB0232440')
        self.assertTrue('description' in attr)
        self.assertEqual(attr['description'],'Contig generated from contig finding genome version 2.5')
        self.assertTrue('Dbxref' in attr)
        self.assertEqual(attr['Dbxref'],'Contig GI Number:90970918,Accession Number:AAFI02000001,SeqID for Genbank:DDB0232440.02')
        self.assertFalse('not_here' in attr)

    def test_percent_decoding(self):
        """Test that percent decoding works on attribute lists
        """
        attributes = "ID=DDB_G0789012;Name=DDB_G0789012_ps;description=putative pseudogene%3B similar to a family of genes%2C including %3Ca href%3D%22%2Fgene%2FDDB_G0234567%22%3EDDB_G0234567%3C%2Fa%3E"
        attr = GFFAttributes(attributes)
        self.assertEqual(attr['ID'],'DDB_G0789012')
        self.assertEqual(attr['Name'],'DDB_G0789012_ps')
        self.assertEqual(attr['description'],'putative pseudogene; similar to a family of genes, including <a href="/gene/DDB_G0234567">DDB_G0234567</a>')

    def test_percent_encoding(self):
        """Test that attributes are properly re-encoded
        """
        dbxref = "Dbxref=Contig GI Number:90970918,Accession Number:AAFI02000001,SeqID for Genbank:DDB0232440.02"
        self.assertEqual(dbxref,str(GFFAttributes(dbxref)))
        parent = "Parent=AF2312,AB2812,abc-3"
        self.assertEqual(parent,str(GFFAttributes(parent)))
        description = "description=putative pseudogene%3B similar to a family of genes%2C including %3Ca href%3D%22%2Fgene%2FDDB_G0234567%22%3EDDB_G0234567%3C%2Fa%3E"
        self.assertEqual(description,str(GFFAttributes(description)))

    def test_recover_representation(self):
        """Test that __repr__ returns original string
        """
        attributes = "ID=DDB0232440;Parent=DDB0232428;Name=DDB0232440;description=Contig generated from contig finding genome version 2.5;Dbxref=Contig GI Number:90970918,Accession Number:AAFI02000001,SeqID for Genbank:DDB0232440.02"
        attr = GFFAttributes(attributes)
        self.assertEqual(attributes,str(attr))

    def test_recover_representation_with_percent_encoding(self):
        """Test that __repr__ returns original string with percent encoding
        """
        attributes = "ID=DDB_G0789012;Name=DDB_G0789012_ps;description=putative pseudogene%3B similar to a family of genes%2C including %3Ca href%3D%22%2Fgene%2FDDB_G0234567%22%3EDDB_G0234567%3C%2Fa%3E"
        attr = GFFAttributes(attributes)
        self.assertEqual(attributes,str(attr))

    def test_recover_representation_with_trailing_semicolon(self):
        """Test that __repr__ returns original string with trailing semicolon
        """
        attributes = "ID=DDB0232440;Parent=DDB0232428;Name=DDB0232440;description=Contig generated from contig finding genome version 2.5;Dbxref=Contig GI Number:90970918,Accession Number:AAFI02000001,SeqID for Genbank:DDB0232440.02;"
        attr = GFFAttributes(attributes)
        self.assertEqual(attributes,str(attr))

    def test_switch_off_encoding(self):
        """Test that __repr__ returns unescaped strings when encoding is off
        """      
        attributes = "ID=DDB_G0789012;Name=DDB_G0789012_ps;description=putative pseudogene%3B similar to a family of genes%2C including %3Ca href%3D%22%2Fgene%2FDDB_G0234567%22%3EDDB_G0234567%3C%2Fa%3E"
        
        attr = GFFAttributes(attributes)
        self.assertTrue(attr.encode())
        self.assertEqual(attributes,str(attr))
        self.assertTrue(attr.encode(True))
        self.assertEqual(attributes,str(attr))
        self.assertFalse(attr.encode(False))
        self.assertEqual("ID=DDB_G0789012;Name=DDB_G0789012_ps;description=putative pseudogene; similar to a family of genes, including <a href=\"/gene/DDB_G0234567\">DDB_G0234567</a>",str(attr))

    def test_empty_attributes_string(self):
        """Test that empty input generates empty output
        """
        attr = GFFAttributes('')
        self.assertEqual('',str(attr))
        attr['ID'] = 'test'
        self.assertEqual('ID=test',str(attr))

    def test_no_null_attributes_for_trailing_semicolon(self):
        """Test trailing ';' in attributes doesn't result in empty 'no-key' attribute
        """
        attr = GFFAttributes("ID=GOAT_ENSP00000332624;")
        self.assertEqual(attr.keys(),['ID'])
        self.assertEqual(attr['ID'],"GOAT_ENSP00000332624")
        self.assertEqual(attr.nokeys(),[])

    def test_nokeys_attributes(self):
        """Test that 'nokeys' attributes are read
        """
        attributes = "ID=GOAT_ENSP00000332668;Shift=2;589008-589009;589535-589539"
        attr = GFFAttributes(attributes)
        self.assertEqual(attr.keys(),['ID','Shift'])
        self.assertEqual(attr['ID'],"GOAT_ENSP00000332668")
        self.assertEqual(attr['Shift'],"2")
        self.assertEqual(attr.nokeys(),['589008-589009','589535-589539'])
        self.assertEqual(str(attr),attributes)

    def test_manipulate_nokeys_attributes(self):
        """Test that 'nokeys' attributes can be manipulated
        """
        attr = GFFAttributes("ID=GOAT_ENSP00000332668;Shift=2;589008-589009;589535-589539")
        self.assertEqual(attr.nokeys(),['589008-589009','589535-589539'])
        nokeys = attr.nokeys()
        nokeys.append('Test')
        self.assertEqual(attr.nokeys(),['589008-589009','589535-589539','Test'])
        del(nokeys[:])
        self.assertEqual(attr.nokeys(),[])

class TestGFFID(unittest.TestCase):
    """Unit tests for GFFID class
    """

    def test_simple_ID(self):
        """Check that name and index can be extracted from ID without code or index
        """
        gffid = GFFID('XYZ123-A')
        self.assertEqual(gffid.code,'')
        self.assertEqual(gffid.name,'XYZ123-A')
        self.assertEqual(gffid.index,0)
        self.assertEqual(str(gffid),'XYZ123-A')        

    def test_full_ID(self):
        """Check that code, name and index can be extracted from full ID
        """
        gffid = GFFID('CDS:XYZ123-A:2')
        self.assertEqual(gffid.code,'CDS')
        self.assertEqual(gffid.name,'XYZ123-A')
        self.assertEqual(gffid.index,2)
        self.assertEqual(str(gffid),'CDS:XYZ123-A:2')

class TestOrderedDictionary(unittest.TestCase):
    """Unit tests for the OrderedDictionary class
    """

    def test_get_and_set(self):
        """Add and retrieve data
        """
        d = OrderedDictionary()
        self.assertEqual(len(d),0)
        d['hello'] = 'goodbye'
        self.assertEqual(d['hello'],'goodbye')

    def test_insert(self):
        """Insert items
        """
        d = OrderedDictionary()
        d['hello'] = 'goodbye'
        self.assertEqual(d.keys(),['hello'])
        self.assertEqual(len(d),1)
        # Insert at start of list
        d.insert(0,'stanley','fetcher')
        self.assertEqual(d.keys(),['stanley','hello'])
        self.assertEqual(len(d),2)
        self.assertEqual(d['stanley'],'fetcher')
        self.assertEqual(d['hello'],'goodbye')
        # Insert in middle
        d.insert(1,'monty','python')
        self.assertEqual(d.keys(),['stanley','monty','hello'])
        self.assertEqual(len(d),3)
        self.assertEqual(d['stanley'],'fetcher')
        self.assertEqual(d['monty'],'python')
        self.assertEqual(d['hello'],'goodbye')

    def test_keeps_things_in_order(self):
        """Check that items are returned in same order as added
        """
        d = OrderedDictionary()
        d['hello'] = 'goodbye'
        d['stanley'] = 'fletcher'
        d['monty'] = 'python'
        self.assertEqual(d.keys(),['hello','stanley','monty'])

    def test_iteration_over_keys(self):
        """Check iterating over keys
        """
        d = OrderedDictionary()
        d['hello'] = 'goodbye'
        d['stanley'] = 'fletcher'
        d['monty'] = 'python'
        try:
            for i in d:
                pass
        except KeyError:
            self.fail("Iteration over OrderedDictionary failed")
