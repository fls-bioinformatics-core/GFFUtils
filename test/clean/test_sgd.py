#!/usr/bin/env python

import unittest
from io import StringIO
from bcftbx.TabFile import TabFile
from GFFUtils.GFFFile import GFFFile
from GFFUtils.clean.sgd import *

class TestGroupByID(unittest.TestCase):

    def setUp(self):
        # List of genes to group
        self.fp = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=CDS:YEL0W02:1;SGD=YEL0W02
chr1\tTest\tCDS\t32611\t34140\t0\t-\t0\tID=CDS:YEL0W02:2;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=CDS:YEL0W03:1;SGD=YEL0W03
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=CDS:YEL0W04:1;SGD=YEL0W04
chr2\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=CDS:YEL0W05:1;SGD=YEL0W05
chr2\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=CDS:YEL0W05:2;SGD=YEL0W05
chr2\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=CDS:YEL0W05:3;SGD=YEL0W05
chr2\tTest\tCDS\t41402\t41831\t0\t-\t0\tID=CDS:YEL0W06:1;SGD=YEL0W06
""")

    def test_group_by_id(self):
        """
        GroupByID: returns features grouped by ID
        """
        gff = GFFFile('test.gff',self.fp)
        groups = GroupByID(gff)
        self.assertEqual(len(groups),6)
        expected_group_ids = ((u"CDS:YEL0W01:1",),
                              (u"CDS:YEL0W02:1",
                               u"CDS:YEL0W02:2",),
                              (u"CDS:YEL0W03:1",),
                              (u"CDS:YEL0W04:1",),
                              (u"CDS:YEL0W05:1",
                               u"CDS:YEL0W05:2",
                               u"CDS:YEL0W05:3",),
                              (u"CDS:YEL0W06:1",))
        for group,expected_group in zip(groups,expected_group_ids):
            for idx,expected_idx in zip(group,expected_group):
                self.assertEqual(idx['attributes']['ID'],
                                 expected_idx)

class TestGFFGetDuplicateSGDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object to read data in
        self.fp = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=YEL0W01;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t32611\t34140\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=YEL0W04;SGD=YEL0W04
chr1\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=YEL0W05;SGD=YEL0W05
chr1\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=YEL0W05;SGD=YEL0W05
chr2\tTest\tCDS\t40406\t40864\t0\t+\t0\tID=YEL0W05;SGD=YEL0W05
chr2\tTest\tCDS\t41402\t41831\t0\t-\t0\tID=YEL0W06;SGD=YEL0W06
""")

    def test_gff_get_duplicate_sgds(self):
        """
        GFFGetDuplicateSGDs: fetch GFF records with duplicated 'SGD' values
        """
        gff = GFFFile('test.gff',self.fp)
        duplicates = GFFGetDuplicateSGDs(gff)
        # Check duplicates
        self.assertEqual(len(duplicates.keys()),2,
                         "wrong number of duplicate SGDs")
        self.assertTrue('YEL0W02' in duplicates.keys())
        self.assertEqual(len(duplicates['YEL0W02']),2)
        self.assertTrue('YEL0W05' in duplicates.keys())
        self.assertEqual(len(duplicates['YEL0W05']),3)

class TestGFFGetDuplicateSGDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object to read data in
        self.fp = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=YEL0W01;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t32611\t34140\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=YEL0W04;SGD=YEL0W04
chr1\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=YEL0W05;SGD=YEL0W05
chr1\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=YEL0W05;SGD=YEL0W05
chr2\tTest\tCDS\t40406\t40864\t0\t+\t0\tID=YEL0W05;SGD=YEL0W05
chr2\tTest\tCDS\t41402\t41831\t0\t-\t0\tID=YEL0W06;SGD=YEL0W06
""")
    
    def test_get_duplicate_sgds(self):
        """Test identifying and returning lines with duplicate SGDs
        """
        gff = GFFFile('test.gff',self.fp)
        duplicates = GFFGetDuplicateSGDs(gff)
        # Check duplicates
        self.assertEqual(len(duplicates.keys()),2,
                         "wrong number of duplicate SGDs")
        self.assertTrue('YEL0W02' in duplicates.keys())
        self.assertEqual(len(duplicates['YEL0W02']),2)
        self.assertTrue('YEL0W05' in duplicates.keys())
        self.assertEqual(len(duplicates['YEL0W05']),3)

class TestGFFResolveDuplicateSGDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object with data for duplicate
        # resolution test
        #
        # Possible duplicate cases are:
        # - unrelated duplicates in same chromosome (YEL0W01)
        # - unrelated duplicates in different chromosomes (YEL0W2)
        # - grouped duplicated in same chromosome (YEL0W03)
        self.fp = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=CDS:YEL0W02:1;SGD=YEL0W02
chr1\tTest\tCDS\t32611\t34140\t0\t-\t0\tID=CDS:YEL0W02:2;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=CDS:YEL0W03:1;SGD=YEL0W03
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=CDS:YEL0W03:2;SGD=YEL0W03
chr1\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=CDS:YEL0W04:1;SGD=YEL0W04
chr1\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
chr2\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=CDS:YEL0W02:1;SGD=YEL0W02
chr2\tTest\tCDS\t41402\t41831\t0\t+\t0\tID=CDS:YEL0W05:1;SGD=YEL0W05
""")
        # Mapping data to resolve all duplicates
        self.mp_resolve_all = StringIO(
u"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with a mapping gene removed
        self.mp_missing_mapping_gene = StringIO(
u"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t-
""")
        # Mapping data with a mapping gene that doesn't match
        # on strand for one duplicate
        self.mp_missing_matching_mapping_gene = StringIO(
u"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t+
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with a mapping gene that doesn't overlap
        self.mp_mapping_gene_no_overlap = StringIO(
u"""YEL0W03\tchr1\t34525\t37004\t-
YEL0W01\tchr1\t41402\t41831\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with multiple mapping genes with same name
        self.mp_multiple_mapping_genes = StringIO(
u"""YEL0W01\tchr1\t28789\t29049\t-
YEL0W03\tchr1\t34525\t37004\t-
YEL0W01\tchr1\t39195\t39569\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
    
    def test_resolve_duplicate_sgds(self):
        """
        GFFResolveDuplicateSGDs: all duplicates can be resolved
        """
        # Load data
        gff = GFFFile('test.gff',self.fp)
        mapping = TabFile('map.txt',self.mp_resolve_all,
                          column_names=('name','chr','start','end','strand'))
        # Fetch duplicates
        duplicates = GFFGetDuplicateSGDs(gff)
        self.assertEqual(len(duplicates),3,"wrong number of duplicates returned")
        # Resolve
        result = GFFResolveDuplicateSGDs(gff,mapping,duplicates,1000)
        # Check results of resolution (GFF data should be unchanged)
        #
        # Selected for discard
        discard = result["discard"]
        self.assertEqual(len(discard),3,"wrong number of duplicates marked for discard")
        self.assertTrue(discard[0] == gff[0])
        self.assertTrue(discard[1] == gff[1])
        self.assertTrue(discard[2] == gff[2])
        #
        # Resolved SGDs
        resolved_sgds = result["resolved_sgds"]
        self.assertEqual(len(resolved_sgds),3)
        self.assertTrue('YEL0W01' in resolved_sgds)
        self.assertTrue('YEL0W02' in resolved_sgds)
        self.assertTrue('YEL0W03' in resolved_sgds)
        #
        # Unresolved SGDs
        self.assertEqual(len(result["unresolved_sgds"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes_after_filter"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_overlaps"]),0)
        self.assertEqual(len(result["unresolved_sgds_multiple_matches"]),0)
    
    def test_resolve_duplicate_sgds_no_mapping_gene(self):
        """
        GFFResolveDuplicateSGDs: missing mapping gene
        """
        # Load data
        gff = GFFFile('test.gff',self.fp)
        mapping = TabFile('map.txt',self.mp_missing_mapping_gene,
                          column_names=('name','chr','start','end','strand'))
        # Fetch duplicates
        duplicates = GFFGetDuplicateSGDs(gff)
        self.assertEqual(len(duplicates),3,"wrong number of duplicates returned")
        # Resolve
        result = GFFResolveDuplicateSGDs(gff,mapping,duplicates,1000)
        # Check results of resolution (GFF data should be unchanged)
        #
        # Selected for discard
        discard = result["discard"]
        self.assertEqual(len(discard),1,"wrong number of duplicates marked for discard")
        self.assertTrue(discard[0] == gff[0])
        #
        # Resolved SGDs
        resolved_sgds = result["resolved_sgds"]
        self.assertEqual(len(resolved_sgds),2)
        self.assertTrue('YEL0W01' in resolved_sgds)
        self.assertTrue('YEL0W03' in resolved_sgds)
        #
        # Unresolved SGDs
        self.assertEqual(len(result["unresolved_sgds"]),1)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes"]),1)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes_after_filter"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_overlaps"]),0)
        self.assertEqual(len(result["unresolved_sgds_multiple_matches"]),0)
        self.assertTrue('YEL0W02' in result["unresolved_sgds_no_mapping_genes"])
        self.assertTrue('YEL0W02' in result["unresolved_sgds"])

    def test_resolve_duplicate_sgds_no_matching_mapping_gene(self):
        """
        GFFResolveDuplicateSGDs: no matching mapping gene
        """
        # Load data
        gff = GFFFile('test.gff',self.fp)
        mapping = TabFile('map.txt',self.mp_missing_matching_mapping_gene,
                          column_names=('name','chr','start','end','strand'))
        # Fetch duplicates
        duplicates = GFFGetDuplicateSGDs(gff)
        self.assertEqual(len(duplicates),3,"wrong number of duplicates returned")
        # Resolve
        result = GFFResolveDuplicateSGDs(gff,mapping,duplicates,1000)
        # Check results of resolution (GFF data should be unchanged)
        #
        # Selected for discard
        discard = result["discard"]
        self.assertEqual(len(discard),3,"wrong number of duplicates marked for discard")
        self.assertTrue(discard[0] == gff[0])
        self.assertTrue(discard[1] == gff[1])
        self.assertTrue(discard[2] == gff[2])
        #
        # Resolved SGDs
        resolved_sgds = result["resolved_sgds"]
        self.assertEqual(len(resolved_sgds),2)
        self.assertTrue('YEL0W01' in resolved_sgds)
        self.assertTrue('YEL0W02' in resolved_sgds)
        #
        # Unresolved SGDs
        self.assertEqual(len(result["unresolved_sgds"]),1)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes_after_filter"]),1)
        self.assertEqual(len(result["unresolved_sgds_no_overlaps"]),0)
        self.assertEqual(len(result["unresolved_sgds_multiple_matches"]),0)
        self.assertTrue('YEL0W03' in result["unresolved_sgds_no_mapping_genes_after_filter"])
        self.assertTrue('YEL0W03' in result["unresolved_sgds"])

    def test_resolve_duplicate_sgds_mapping_gene_with_no_overlap(self):
        """
        GFFResolveDuplicateSGDs: mapping gene with no overlap
        """
        # Load data
        gff = GFFFile('test.gff',self.fp)
        mapping = TabFile('map.txt',self.mp_mapping_gene_no_overlap,
                          column_names=('name','chr','start','end','strand'))
        # Fetch duplicates
        duplicates = GFFGetDuplicateSGDs(gff)
        self.assertEqual(len(duplicates),3,"wrong number of duplicates returned")
        # Resolve
        result = GFFResolveDuplicateSGDs(gff,mapping,duplicates,1000)
        # Check results of resolution (GFF data should be unchanged)
        #
        # Selected for discard
        discard = result["discard"]
        self.assertEqual(len(discard),2,"wrong number of duplicates marked for discard")
        self.assertTrue(discard[0] == gff[1])
        self.assertTrue(discard[1] == gff[2])
        #
        # Resolved SGDs
        resolved_sgds = result["resolved_sgds"]
        self.assertEqual(len(resolved_sgds),2)
        self.assertTrue('YEL0W02' in resolved_sgds)
        self.assertTrue('YEL0W03' in resolved_sgds)
        #
        # Unresolved SGDs
        self.assertEqual(len(result["unresolved_sgds"]),1)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes_after_filter"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_overlaps"]),1)
        self.assertEqual(len(result["unresolved_sgds_multiple_matches"]),0)
        self.assertTrue('YEL0W01' in result["unresolved_sgds_no_overlaps"])
        self.assertTrue('YEL0W01' in result["unresolved_sgds"])

    def test_resolve_duplicate_sgds_multiple_mapping_genes_with_same_name(self):
        """
        GFFResolveDuplicateSGDs: multiple mapping genes with the same name
        """
        # Load data
        gff = GFFFile('test.gff',self.fp)
        mapping = TabFile('map.txt',self.mp_multiple_mapping_genes,
                          column_names=('name','chr','start','end','strand'))
        # Fetch duplicates
        duplicates = GFFGetDuplicateSGDs(gff)
        self.assertEqual(len(duplicates),3,"wrong number of duplicates returned")
        # Resolve
        result = GFFResolveDuplicateSGDs(gff,mapping,duplicates,1000)
        # Check results of resolution (GFF data should be unchanged)
        #
        # Selected for discard
        discard = result["discard"]
        self.assertEqual(len(discard),2,"wrong number of duplicates marked for discard")
        self.assertTrue(discard[0] == gff[1])
        self.assertTrue(discard[1] == gff[2])
        #
        # Resolved SGDs
        resolved_sgds = result["resolved_sgds"]
        self.assertEqual(len(resolved_sgds),2)
        self.assertTrue('YEL0W02' in resolved_sgds)
        self.assertTrue('YEL0W03' in resolved_sgds)
        #
        # Unresolved SGDs
        self.assertEqual(len(result["unresolved_sgds"]),1)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_mapping_genes_after_filter"]),0)
        self.assertEqual(len(result["unresolved_sgds_no_overlaps"]),0)
        self.assertEqual(len(result["unresolved_sgds_multiple_matches"]),1)
        self.assertTrue('YEL0W01' in result["unresolved_sgds_multiple_matches"])
        self.assertTrue('YEL0W01' in result["unresolved_sgds"])

class TestGFFGroupSGDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object to read data in
        self.fp = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=YEL0W01;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t32611\t34140\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=YEL0W04;SGD=YEL0W04
chr1\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=YEL0W01;SGD=YEL0W01
""")
        # The expected ID assignments for the input above
        self.expected_ids = ("CDS:YEL0W01:1",
                             "CDS:YEL0W02:1",
                             "CDS:YEL0W02:2",
                             "CDS:YEL0W03:1",
                             "CDS:YEL0W03:2",
                             "CDS:YEL0W04:1",
                             "CDS:YEL0W03:3",
                             "CDS:YEL0W01:1")

    def test_gff_group_sgds(self):
        """
        GFFGroupSGDs: relabel ID with SGD information
        """
        gff = GFFFile('test.gff',self.fp)
        GFFGroupSGDs(gff)
        for line,expected_id in zip(gff,self.expected_ids):
            self.assertEqual(line['attributes']['ID'],expected_id)

class TestGFFInsertMissingGenes(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF data
        self.fp = StringIO(
u"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=CDS:YEL0W02:1;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=CDS:YEL0W04:1;SGD=YEL0W04
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=CDS:YEL0W04:2;SGD=YEL0W04
chr2\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=CDS:YEL0W05:1;SGD=YEL0W05
chr2\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=CDS:YEL0W06:1;SGD=YEL0W06
chr2\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=CDS:YEL0W06:2;SGD=YEL0W06
""")
        # Make a file-like object for mapping data
        self.mp = StringIO(
u"""YEL0W03\tchr1\t32611\t34140\t-
YEL0W06\tchr2\t49195\t49569\t-
""")

    def test_gff_insert_missing_genes(self):
        """
        GFFInsertMissingGenes: insert gene from mapping data file
        """
        gff = GFFFile('test.gff',self.fp)
        mapping = TabFile('map.txt',self.mp,
                          column_names=('name','chr','start','end','strand'))
        # Insert missing genes
        GFFInsertMissingGenes(gff,mapping)
        # Check: YEL0W03 should be inserted at index 2
        i = 2
        idx = gff[i]['attributes']['ID']
        self.assertEqual(idx,'CDS:YEL0W03:1',"incorrect ID at position %d" % i)
        # Check: YEL0W06 should not have been inserted
        for i in range(len(gff)):
            idx = gff[i]['attributes']['ID']
            self.assertNotEqual(idx,'CDS:YEL0W06:3',"wrong ID at position %d" %i)
        # Check: no leading ';' on the string representation
        self.assertNotEqual(str(gff[i]['attributes'])[0],';',"Erroneous leading semicolon: %s"
                            % str(gff[i]['attributes']))
