#!/usr/bin/env python

import unittest
from io import StringIO
from GFFUtils.GFFFile import GFFFile
from GFFUtils.clean.generic import *

class TestGFFUpdateAttributes(unittest.TestCase):

    def setUp(self):
        # Make file-like object to read data in
        self.fp = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=abc;kaks=-le+100;SGD=YEL0W;ncbi=-1e+100;Name=def;
""")
        self.fp2 = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=abc;kaks=-le+100;SGD=YEL0W;ncbi=-1e+100;Name=def;123-234;456-567;
""")
        self.fp3 = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=abc;kaks=-le+100;SGD=;ncbi=-1e+100;Name=def;
""")

    def test_update_attributes_exclude_keys(self):
        """
        GFFUpdateAttributes: exclude specific attributes
        """
        gff = GFFFile('test.gff',self.fp)
        # Check that attributes are present initially
        attributes = gff[0]['attributes']
        for attr in ['ID','Name','SGD','ncbi','kaks']:
            self.assertTrue(attr in attributes.keys())
        # Do the exclusion operation
        GFFUpdateAttributes(gff,exclude_keys=['ncbi','kaks'])
        # Check that expected attributes have been removed
        attributes = gff[0]['attributes']
        for attr in ['ID','Name','SGD']:
            self.assertTrue(attr in attributes.keys())
        for attr in ['ncbi','kaks']:
            self.assertTrue(attr not in attributes.keys())

    def test_update_attributes_exclude_nokeys(self):
        """
        GFFUpdateAttributes: exclude 'nokeys' attributes
        """
        gff = GFFFile('test.gff',self.fp2)
        # Check that nokey attributes are present initially
        attributes = gff[0]['attributes']
        self.assertEqual(attributes.nokeys(),['123-234','456-567'])
        # Do the exclusion operation
        GFFUpdateAttributes(gff,exclude_nokeys=True)
        # Check that nokey attributes have been removed
        attributes = gff[0]['attributes']
        self.assertEqual(attributes.nokeys(),[])

    def test_update_attributes_replace_values(self):
        """
        GFFUpdateAttributes: replace specific attributes
        """
        gff = GFFFile('test.gff',self.fp)
        GFFUpdateAttributes(gff,update_keys={'ID':'SGD','Name':'SGD'})
        # Check that attributes have the expected values
        attributes = gff[0]['attributes']
        sgd = attributes['SGD']
        for attr in ['ID','Name']:
            self.assertTrue(attributes[attr] == sgd)

    def test_gff_update_attributes_replace_with_empty_value(self):
        """
        GFFUpdateAttributes: replace attribute with empty value
        """
        gff = GFFFile("test.gff",fp=self.fp3)
        GFFUpdateAttributes(gff,update_keys={ 'Name':'SGD' },
                            no_empty_values=False)
        # Check that attributes have the expected values
        attributes = gff[0]['attributes']
        self.assertEqual(attributes['Name'],'')

    def test_gff_update_attributes_update_dont_allow_empty_value(self):
        """
        GFFUpdateAttributes: don't allow attribute replacement with empty value
        """
        gff = GFFFile("test.gff",fp=self.fp3)
        GFFUpdateAttributes(gff,update_keys={ 'Name':'SGD' },
                            no_empty_values=True)
        # Check that attributes have the expected values
        attributes = gff[0]['attributes']
        self.assertEqual(attributes['Name'],'def')

class TestGFFAddExonIDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF pseudo-data
        self.fp = StringIO(
u"""chr1\tTest\texon\t1890\t3287\t.\t+\t.\tParent=DDB0216437
chr1\tTest\texon\t3848\t4855\t.\t+\t.\tParent=DDB0216438
chr1\tTest\tCDS\t5505\t7769\t.\t+\t.\tParent=DDB0216439
chr1\tTest\tCDS\t8308\t9522\t.\t-\t.\tParent=DDB0216440
chr1\tTest\texon\t9635\t9889\t.\t-\t.\tParent=DDB0216441
chr1\tTest\texon\t10033\t11199\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t11264\t11952\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t12069\t12183\t.\t+\t.\tParent=DDB0216442
chr1\tTest\tgene\t12436\t13044\t.\t-\t.\tParent=DDB0216443
chr1\tTest\texon\t17379\t17386\t.\t+\t.\tParent=DDB0216445
""")

    def test_gff_add_exon_ids(self):
        """
        GFFAddExonIDs: constructs and adds missing ID attributes for exons
        """
        gff = GFFFile('test.gff',self.fp)
        # Add the exon IDs
        gff = GFFAddExonIDs(gff)
        # Check that all exons have an ID
        for data in gff:
            attr = data['attributes']
            if data['feature'] == 'exon':
                self.assertTrue('ID' in attr,"No ID attribute found")
                self.assertEqual(attr.keys()[0],'ID',
                                 "ID attribute should be first")
            else:
                self.assertFalse('ID' in attr,"ID shouldn't be present")

class TestGFFAddIDAttributes(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF pseudo-data
        self.fp = StringIO(
u"""chr1\tTest\texon\t1890\t3287\t.\t+\t.\tParent=DDB0216437
chr1\tTest\texon\t3848\t4855\t.\t+\t.\tParent=DDB0216438
chr1\tTest\tCDS\t5505\t7769\t.\t+\t.\tParent=DDB0216439
chr1\tTest\tCDS\t8308\t9522\t.\t-\t.\tParent=DDB0216440
chr1\tTest\texon\t9635\t9889\t.\t-\t.\tParent=DDB0216441
chr1\tTest\texon\t10033\t11199\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t11264\t11952\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t12069\t12183\t.\t+\t.\tParent=DDB0216442
chr1\tTest\tgene\t12436\t13044\t.\t-\t.\tID=DDB012345678;Parent=DDB0216443
chr1\tTest\texon\t17379\t17386\t.\t+\t.\tParent=DDB0216445
""")
        self.expected_ids = ("exon:DDB0216437:00000001",
                             "exon:DDB0216438:00000002",
                             "CDS:DDB0216439:00000003",
                             "CDS:DDB0216440:00000004",
                             "exon:DDB0216441:00000005",
                             "exon:DDB0216442:00000006",
                             "exon:DDB0216442:00000007",
                             "exon:DDB0216442:00000008",
                             "DDB012345678",
                             "exon:DDB0216445:00000009",)

    def test_gff_add_id_attributes(self):
        """
        GFFAddIDAttributes: constructs and adds missing ID attributes
        """
        gff = GFFFile('test.gff',self.fp)
        # Add the exon IDs
        gff = GFFAddIDAttributes(gff)
        # Check that all features have an ID
        for data,expected_id in zip(gff,self.expected_ids):
            attr = data['attributes']
            self.assertTrue('ID' in attr,"No ID attribute found")
            self.assertEqual(attr.keys()[0],'ID',
                             "ID attribute should be first")
            self.assertEqual(attr['ID'],expected_id)

class TestGFFDecodeAttributes(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF pseudo-data
        self.fp = StringIO(
u"""chr1\t.\tgene\t5505\t7769\t.\t+\t.\tID=DDB_G0267182;Name=DDB_G0123456;description=ORF2 protein fragment of DIRS1 retrotransposon%3B refer to Genbank M11339 for full-length element
chr1\t.\tgene\t21490\t23468\t.\t+\t.\tID=DDB_G0267204;Name=DDB_G0123456;description=putative pseudogene%3B similar to a family of genes%2C including %3Ca href%3D%22%2Fgene%2FDDB_G0267252%22%3EDDB_G0267252%3C%2Fa%3E
""")

    def test_gff_decode_attributes(self):
        """
        GFFDecodeAttributes: re-encode escape sequences
        """
        gff = GFFFile('test.gff',self.fp)
        # Decode the attribute data
        gff = GFFDecodeAttributes(gff)
        # Check decoding
        self.assertEqual("%s" % gff[0]['attributes'],"ID=DDB_G0267182;Name=DDB_G0123456;description=ORF2 protein fragment of DIRS1 retrotransposon; refer to Genbank M11339 for full-length element")
        self.assertEqual("%s" % gff[1]['attributes'],"ID=DDB_G0267204;Name=DDB_G0123456;description=putative pseudogene; similar to a family of genes, including <a href=\"/gene/DDB_G0267252\">DDB_G0267252</a>")


