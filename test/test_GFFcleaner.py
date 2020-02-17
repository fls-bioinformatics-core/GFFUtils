#!/usr/bin/env python

import unittest
from io import StringIO
from GFFUtils.GFFcleaner import *

class TestGroupGeneSubsets(unittest.TestCase):

    def setUp(self):
        # List of genes to group
        self.fp = StringIO(
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=CDS:YEL0W02:1;SGD=YEL0W02
chr1\tTest\tCDS\t32611\t34140\t0\t-\t0\tID=CDS:YEL0W02:2;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=CDS:YEL0W03:1;SGD=YEL0W03
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=CDS:YEL0W04:1;SGD=YEL0W04
chr2\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=CDS:YEL0W05:1;SGD=YEL0W05
chr2\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=CDS:YEL0W05:2;SGD=YEL0W05
chr2\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=CDS:YEL0W05:3;SGD=YEL0W05
chr2\tTest\tCDS\t41402\t41831\t0\t-\t0\tID=CDS:YEL0W06:1;SGD=YEL0W06
""")

    def test_group_gene_subsets_by_SGD(self):
        """Test grouping of genes by SGD attribute
        """
        gff = GFFFile('test.gff',self.fp)
        groups = GroupGeneSubsets(gff)
        self.assertEqual(len(groups),6,"expected 6 groups, got %d" % len(groups))

class TestGFFUpdateAttributes(unittest.TestCase):

    def setUp(self):
        # Make file-like object to read data in
        self.fp = StringIO(
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=abc;kaks=-le+100;SGD=YEL0W;ncbi=-1e+100;Name=def;
""")
        self.fp2 = StringIO(
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=abc;kaks=-le+100;SGD=YEL0W;ncbi=-1e+100;Name=def;123-234;456-567;
""")

    def test_update_attributes_exclude_keys(self):
        """Test excluding specific attributes from the attribute field
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
        """Test excluding 'nokeys' attributes from the attribute field
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
        """Test replacing specific attributes from the attribute field
        """
        gff = GFFFile('test.gff',self.fp)
        GFFUpdateAttributes(gff,update_keys={'ID':'SGD','Name':'SGD'})
        # Check that attributes have the expected values
        attributes = gff[0]['attributes']
        sgd = attributes['SGD']
        for attr in ['ID','Name']:
            self.assertTrue(attributes[attr] == sgd)

class TestGFFGetDuplicateSGDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object to read data in
        self.fp = StringIO(
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=YEL0W01;SGD=YEL0W01
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
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
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
"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with a mapping gene removed
        self.mp_missing_mapping_gene = StringIO(
"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t-
""")
        # Mapping data with a mapping gene that doesn't match
        # on strand for one duplicate
        self.mp_missing_matching_mapping_gene = StringIO(
"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t+
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with a mapping gene that doesn't overlap
        self.mp_mapping_gene_no_overlap = StringIO(
"""YEL0W03\tchr1\t34525\t37004\t-
YEL0W01\tchr1\t41402\t41831\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with multiple mapping genes with same name
        self.mp_multiple_mapping_genes = StringIO(
"""YEL0W01\tchr1\t28789\t29049\t-
YEL0W03\tchr1\t34525\t37004\t-
YEL0W01\tchr1\t39195\t39569\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
    
    def test_resolve_duplicate_sgds(self):
        """Test resolving duplicate SGDs (all duplicates can be resolved)
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
        """Test resolving duplicate SGDs (missing mapping gene)
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
        """Test resolving duplicate SGDs (no matching mapping gene)
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
        """Test resolving duplicate SGDs (mapping gene with no overlap)
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
        """Test resolving duplicate SGDs (multiple mapping genes with the same name)
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
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=YEL0W01;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t32611\t34140\t0\t-\t0\tID=YEL0W02;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=YEL0W04;SGD=YEL0W04
chr1\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=YEL0W03;SGD=YEL0W03
chr1\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=YEL0W01;SGD=YEL0W01
""")
        # The expected ID assignments for the input above
        self.ids = ["CDS:YEL0W01:1",
                    "CDS:YEL0W02:1",
                    "CDS:YEL0W02:2",
                    "CDS:YEL0W03:1",
                    "CDS:YEL0W03:2",
                    "CDS:YEL0W04:1",
                    "CDS:YEL0W03:3",
                    "CDS:YEL0W01:1"]

    def test_gff_group_sgds(self):
        """Test ID attributes are correctly assigned
        """
        gff = GFFFile('test.gff',self.fp)
        # Group by SGD
        GFFGroupSGDs(gff)
        # Check the ID attribute for each line
        for i in range(len(gff)):
            idx = gff[i]['attributes']['ID']
            self.assertEqual(idx,self.ids[i],"incorrect ID at position %d" % i)

class TestGFFInsertMissingGenes(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF data
        self.fp = StringIO(
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=CDS:YEL0W02:1;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=CDS:YEL0W04:1;SGD=YEL0W04
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=CDS:YEL0W04:2;SGD=YEL0W04
chr2\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=CDS:YEL0W05:1;SGD=YEL0W05
chr2\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=CDS:YEL0W06:1;SGD=YEL0W06
chr2\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=CDS:YEL0W06:2;SGD=YEL0W06
""")
        # Make a file-like object for mapping data
        self.mp = StringIO(
"""YEL0W03\tchr1\t32611\t34140\t-
YEL0W06\tchr2\t49195\t49569\t-
"""
)

    def test_insert_missing_gene(self):
        """Test inserting a missing gene into GFF data
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

class TestGFFAddExonIDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF pseudo-data
        self.fp = StringIO(
"""chr1\tTest\texon\t1890\t3287\t.\t+\t.\tParent=DDB0216437
chr1\tTest\texon\t3848\t4855\t.\t+\t.\tParent=DDB0216438
chr1\tTest\tCDS\t5505\t7769\t.\t+\t.\tParent=DDB0216439
chr1\tTest\tCDS\t8308\t9522\t.\t-\t.\tParent=DDB0216440
chr1\tTest\texon\t9635\t9889\t.\t-\t.\tParent=DDB0216441
chr1\tTest\texon\t10033\t11199\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t11264\t11952\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t12069\t12183\t.\t+\t.\tParent=DDB0216442
chr1\tTest\tgene\t12436\t13044\t.\t-\t.\tParent=DDB0216443
chr1\tTest\texon\t17379\t17386\t.\t+\t.	Parent=DDB0216445
""")

    def test_add_exon_id(self):
        """Test adding ID attributes to exon records
        """
        gff = GFFFile('test.gff',self.fp)
        # Add the exon IDs
        gff = GFFAddExonIDs(gff)
        # Check that all exons have an ID
        for data in gff:
            if data['feature'] == 'exon':
                attr = data['attributes']
                self.assertTrue('ID' in attr,"No ID attribute found")
                self.assertEqual(attr.keys()[0],'ID',"ID attribute should be first")

class TestGFFAddIDAttributes(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF pseudo-data
        self.fp = StringIO(
"""chr1\tTest\texon\t1890\t3287\t.\t+\t.\tParent=DDB0216437
chr1\tTest\texon\t3848\t4855\t.\t+\t.\tParent=DDB0216438
chr1\tTest\tCDS\t5505\t7769\t.\t+\t.\tParent=DDB0216439
chr1\tTest\tCDS\t8308\t9522\t.\t-\t.\tParent=DDB0216440
chr1\tTest\texon\t9635\t9889\t.\t-\t.\tParent=DDB0216441
chr1\tTest\texon\t10033\t11199\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t11264\t11952\t.\t+\t.\tParent=DDB0216442
chr1\tTest\texon\t12069\t12183\t.\t+\t.\tParent=DDB0216442
chr1\tTest\tgene\t12436\t13044\t.\t-\t.\tID=DDB012345678;Parent=DDB0216443
chr1\tTest\texon\t17379\t17386\t.\t+\t.	Parent=DDB0216445
""")

    def test_add_ids(self):
        """Test adding ID attributes
        """
        gff = GFFFile('test.gff',self.fp)
        # Add the exon IDs
        gff = GFFAddIDAttributes(gff)
        # Check that all exons have an ID
        for data in gff:
            attr = data['attributes']
            self.assertTrue('ID' in attr,"No ID attribute found")
            self.assertEqual(attr.keys()[0],'ID',"ID attribute should be first")

class TestGFFDecodeAttributes(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF pseudo-data
        self.fp = StringIO(
"""chr1\t.\tgene\t5505\t7769\t.\t+\t.\tID=DDB_G0267182;Name=DDB_G0123456;description=ORF2 protein fragment of DIRS1 retrotransposon%3B refer to Genbank M11339 for full-length element
chr1\t.\tgene\t21490\t23468\t.\t+\t.\tID=DDB_G0267204;Name=DDB_G0123456;description=putative pseudogene%3B similar to a family of genes%2C including %3Ca href%3D%22%2Fgene%2FDDB_G0267252%22%3EDDB_G0267252%3C%2Fa%3E
""")

    def test_decode_attributes(self):
        """Test adding ID attributes
        """
        gff = GFFFile('test.gff',self.fp)
        # Decode the attribute data
        gff = GFFDecodeAttributes(gff)
        # Check decoding
        self.assertEqual("%s" % gff[0]['attributes'],"ID=DDB_G0267182;Name=DDB_G0123456;description=ORF2 protein fragment of DIRS1 retrotransposon; refer to Genbank M11339 for full-length element")
        self.assertEqual("%s" % gff[1]['attributes'],"ID=DDB_G0267204;Name=DDB_G0123456;description=putative pseudogene; similar to a family of genes, including <a href=\"/gene/DDB_G0267252\">DDB_G0267252</a>")
