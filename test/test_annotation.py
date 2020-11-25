#!/usr/bin/env python

import unittest
import tempfile
import shutil
import os
from io import StringIO
from GFFUtils.GFFFile import GFFFile
from GFFUtils.GFFFile import GFFDataLine
from GFFUtils.GTFFile import GTFFile
from GFFUtils.GTFFile import GTFDataLine
from GFFUtils.annotation import *

# Example GFF file fragment
gff_data = u"""##gff-version   3
##sequence-region DDB0232428 1 4923596
##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=352472
# generated: Wed Jan 25 03:46:58 2012
DDB0232429	Sequencing Center	mRNA	6679320	6680012	.	+	.	ID=DDB0166998;Parent=DDB_G0276345;Name=DDB0166998;
DDB0232429	Sequencing Center	exon	6679320	6679397	.	+	.	Parent=DDB0166998
DDB0232429	Sequencing Center	CDS	6679320	6679397	.	+	.	Parent=DDB0166998
DDB0232429	Sequencing Center	exon	6679470	6679535	.	+	.	Parent=DDB0166998
DDB0232429	Sequencing Center	CDS	6679470	6679535	.	+	.	Parent=DDB0166998
DDB0232429	Sequencing Center	exon	6679617	6680012	.	+	.	Parent=DDB0166998
DDB0232429	Sequencing Center	CDS	6679617	6680012	.	+	.	Parent=DDB0166998
DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=Description of gene naa20
DDB0232429	Sequencing Center	mRNA	6679320	6680012	.	+	.	ID=DDB0166998;Parent=DDB_G0276345;Name=DDB0166998;
DDB0232429	dictyBase Curator	mRNA	6679320	6680012	.	+	.	ID=DDB0238097;Parent=DDB_G0276345;Name=DDB0238097;
DDB0232429	Sequencing Center	mRNA	5954835	5955486	.	+	.	ID=DDB0167147;Parent=DDB_G0275629;Name=DDB0167147;
DDB0232429	Sequencing Center	exon	5954835	5954850	.	+	.	Parent=DDB0167147
DDB0232429	Sequencing Center	CDS	5954835	5954850	.	+	.	Parent=DDB0167147
DDB0232429	Sequencing Center	exon	5955332	5955486	.	+	.	Parent=DDB0167147
DDB0232429	Sequencing Center	CDS	5955332	5955486	.	+	.	Parent=DDB0167147
DDB0232429	.	gene	5954835	5955486	.	+	.	ID=DDB_G0275629;Name=DDB_G0275629;description=Description of gene DDB_G0275629
"""

# Example GTF fragment
gtf_data = u"""##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)
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
"""

class TestGFFAnnotationLookup(unittest.TestCase):

    def test_gff_annotation_lookup_from_gff(self):
        """
        GFFAnnotationLookup: lookup from GFF file
        """
        # Load GFF data
        gff = GFFFile("test.gff",StringIO(gff_data))
        lookup = GFFAnnotationLookup(gff)
        # getDataFromID
        data = lookup.getDataFromID("DDB_G0276345")
        self.assertEqual(len(data),1)
        self.assertEqual(str(data[0]),
                         str(GFFDataLine("DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=Description of gene naa20")))
        # getParentID
        self.assertEqual(lookup.getParentID("DDB0166998"),"DDB_G0276345")
        # getAncestorGene
        self.assertEqual(str(lookup.getAncestorGene("DDB0166998")),
                         "DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=Description of gene naa20")
        # getAnnotation
        annot = lookup.getAnnotation("DDB0166998")
        self.assertEqual(annot.parent_feature_name,"DDB0166998")
        self.assertEqual(annot.parent_feature_type,"mRNA")
        self.assertEqual(annot.parent_feature_parent,"DDB_G0276345")
        self.assertEqual(annot.parent_gene_name,"naa20")
        self.assertEqual(annot.chr,"DDB0232429")
        self.assertEqual(annot.start,6679320)
        self.assertEqual(annot.end,6680012)
        self.assertEqual(annot.strand,"+")
        self.assertEqual(annot.gene_locus,"DDB0232429:6679320-6680012")
        self.assertEqual(annot.gene_length,692)
        self.assertEqual(annot.description,"Description of gene naa20")

    def test_gff_annotation_lookup_from_gtf(self):
        """
        GFFAnnotationLookup: lookup from GTF file
        """
        # Load GFF data
        gtf = GTFFile("test.gtf",StringIO(gtf_data))
        lookup = GFFAnnotationLookup(gtf)
        # getDataFromID
        data = lookup.getDataFromID("ENSG00000223972.4")
        self.assertEqual(len(data),1)
        self.assertEqual(str(data[0]),
                         str(GTFDataLine("""chr1	HAVANA	gene	11869	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";""")))
        # getParentID
        self.assertRaises(KeyError,
                          lookup.getParentID,
                          "ENSG00000223972.4")
        # getAncestorGene
        self.assertEqual(lookup.getAncestorGene("ENSG00000223972.4"),None)
        # getAnnotation
        annot = lookup.getAnnotation("ENSG00000223972.4")
        self.assertEqual(annot.parent_feature_name,"ENSG00000223972.4")
        self.assertEqual(annot.parent_feature_type,"gene")
        self.assertEqual(annot.parent_feature_parent,None)
        self.assertEqual(annot.parent_gene_name,"DDX11L1")
        self.assertEqual(annot.chr,"chr1")
        self.assertEqual(annot.start,11869)
        self.assertEqual(annot.end,14412)
        self.assertEqual(annot.strand,"+")
        self.assertEqual(annot.gene_locus,"chr1:11869-14412")
        self.assertEqual(annot.gene_length,2543)
        self.assertEqual(annot.description,"")

class TestGFFAnnotation(unittest.TestCase):

    def test_empty_gff_annotation(self):
        """
        GFFAnnotation: no data
        """
        annot = GFFAnnotation("DDB0166998")
        self.assertEqual(annot.parent_feature_name,"DDB0166998")
        self.assertEqual(annot.parent_feature_type,"")
        self.assertEqual(annot.parent_feature_parent,"")
        self.assertEqual(annot.parent_gene_name,"")
        self.assertEqual(annot.gene_locus,"")
        self.assertEqual(annot.description,"")
        self.assertEqual(annot.chr,"")
        self.assertEqual(annot.start,"")
        self.assertEqual(annot.end,"")
        self.assertEqual(annot.strand,"")
        self.assertEqual(annot.gene_length,"")

    def test_gff_annotation_set_manually(self):
        """
        GFFAnnotation: data set manually
        """
        annot = GFFAnnotation("DDB0166998")
        annot.parent_feature_type = "mRNA"
        annot.parent_feature_parent = "DDB_G0276345"
        annot.parent_gene_name = "naa20"
        annot.chr = "DDB0232429"
        annot.start = 6679320
        annot.end = 6680012
        annot.strand = "+"
        annot.description = "Description of gene naa20"
        self.assertEqual(annot.parent_feature_name,"DDB0166998")
        self.assertEqual(annot.parent_feature_type,"mRNA")
        self.assertEqual(annot.parent_feature_parent,"DDB_G0276345")
        self.assertEqual(annot.parent_gene_name,"naa20")
        self.assertEqual(annot.gene_locus,"DDB0232429:6679320-6680012")
        self.assertEqual(annot.description,"Description of gene naa20")
        self.assertEqual(annot.chr,"DDB0232429")
        self.assertEqual(annot.start,6679320)
        self.assertEqual(annot.end,6680012)
        self.assertEqual(annot.strand,"+")
        self.assertEqual(annot.gene_length,692)

    def test_gff_annotation_set_from_feature(self):
        """
        GFFAnnotation: data set from feature
        """
        feature = GFFDataLine("DDB0232429	Sequencing Center	mRNA	6679320	6680012	.	+	.	ID=DDB0166998;Parent=DDB_G0276345;Name=DDB0166998;")
        annot = GFFAnnotation("DDB0166998",feature)
        annot.parent_gene_name = "naa20"
        annot.chr = "DDB0232429"
        annot.start = 6679320
        annot.end = 6680012
        annot.strand = "+"
        annot.description = "Description of gene naa20"
        self.assertEqual(annot.parent_feature_name,"DDB0166998")
        self.assertEqual(annot.parent_feature_type,"mRNA")
        self.assertEqual(annot.parent_feature_parent,"DDB_G0276345")
        self.assertEqual(annot.parent_gene_name,"naa20")
        self.assertEqual(annot.gene_locus,"DDB0232429:6679320-6680012")
        self.assertEqual(annot.description,"Description of gene naa20")
        self.assertEqual(annot.chr,"DDB0232429")
        self.assertEqual(annot.start,6679320)
        self.assertEqual(annot.end,6680012)
        self.assertEqual(annot.strand,"+")
        self.assertEqual(annot.gene_length,692)

    def test_gff_annotation_set_from_feature_and_gene(self):
        """
        GFFAnnotation: data set from feature and associated gene
        """
        feature = GFFDataLine("DDB0232429	Sequencing Center	mRNA	6679320	6680012	.	+	.	ID=DDB0166998;Parent=DDB_G0276345;Name=DDB0166998;")
        gene = GFFDataLine("DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=Description of gene naa20")
        annot = GFFAnnotation("DDB0166998",feature,gene)
        self.assertEqual(annot.parent_feature_name,"DDB0166998")
        self.assertEqual(annot.parent_feature_type,"mRNA")
        self.assertEqual(annot.parent_feature_parent,"DDB_G0276345")
        self.assertEqual(annot.parent_gene_name,"naa20")
        self.assertEqual(annot.gene_locus,"DDB0232429:6679320-6680012")
        self.assertEqual(annot.description,"Description of gene naa20")
        self.assertEqual(annot.chr,"DDB0232429")
        self.assertEqual(annot.start,6679320)
        self.assertEqual(annot.end,6680012)
        self.assertEqual(annot.strand,"+")
        self.assertEqual(annot.gene_length,692)

class TestHTSeqCountFile(unittest.TestCase):

    def setUp(self):
        # Temporary directory
        self.wd = tempfile.mkdtemp()
        # htseq-count fragment
        self.htseq_count = """DDB0166998	2
DDB0167147	0
DDB0167277	16
no_feature	120515
ambiguous	422601
too_low_aQual	0
not_aligned	19996043
alignment_not_unique	0
"""
        # Write to temporary file
        self.htseq_count_file = os.path.join(self.wd,
                                             "htseq_counts.txt")
        with open(self.htseq_count_file,'wt') as fp:
            fp.write(self.htseq_count)

    def tearDown(self):
        if os.path.exists(self.wd):
            shutil.rmtree(self.wd)

    def test_htseq_count_file(self):
        """
        HTSeqCountFile: retrieve htseq-count data from file
        """
        htseq = HTSeqCountFile(self.htseq_count_file)
        # feature_IDs
        self.assertEqual(htseq.feature_IDs(),['DDB0166998',
                                              'DDB0167147',
                                              'DDB0167277'])
        # count
        self.assertEqual(htseq.count('DDB0166998'),2)
        self.assertEqual(htseq.count('DDB0167147'),0)
        self.assertEqual(htseq.count('DDB0167277'),16)
        # table
        self.assertEqual(htseq.table().keys(),["total_counted_into_genes",
                                               "no_feature",
                                               "ambiguous",
                                               "too_low_aQual",
                                               "not_aligned",
                                               "alignment_not_unique"])
        self.assertEqual(htseq.table()["total_counted_into_genes"],18)
        self.assertEqual(htseq.table()["no_feature"],120515)
        self.assertEqual(htseq.table()["ambiguous"],422601)
        self.assertEqual(htseq.table()["too_low_aQual"],0)
        self.assertEqual(htseq.table()["not_aligned"],19996043)
        self.assertEqual(htseq.table()["alignment_not_unique"],0)

class TestAnnotateFeatureData(unittest.TestCase):

    def setUp(self):
        # Temporary directory
        self.wd = tempfile.mkdtemp()
        # Feature data fragment
        self.feature_data = """#Gene
DDB0166998
DDB0167147
"""
        # Write to temporary file
        self.feature_data_file = os.path.join(self.wd,
                                             "features.txt")
        with open(self.feature_data_file,'wt') as fp:
            fp.write(self.feature_data)
        # Output file name
        self.out_file = os.path.join(self.wd,"out.txt")

    def tearDown(self):
        if os.path.exists(self.wd):
            shutil.rmtree(self.wd)

    def test_annotate_feature_data(self):
        """
        annotate_feature_data: annotates from GFF file
        """
        # Load GFF data
        gff = GFFFile("test.gff",StringIO(gff_data))
        lookup = GFFAnnotationLookup(gff)
        # Do annotation
        annotate_feature_data(lookup,
                              self.feature_data_file,
                              self.out_file)
        # Check output
        with open(self.out_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             """Gene	exon_parent	feature_type_exon_parent	gene_ID	gene_name	chr	start	end	strand	gene_length	locus	description
DDB0166998	DDB0166998	mRNA	DDB_G0276345	naa20	DDB0232429	6679320	6680012	+	692	DDB0232429:6679320-6680012	Description of gene naa20
DDB0167147	DDB0167147	mRNA	DDB_G0275629	DDB_G0275629	DDB0232429	5954835	5955486	+	651	DDB0232429:5954835-5955486	Description of gene DDB_G0275629
""")

class TestAnnotateHtseqCountData(unittest.TestCase):

    def setUp(self):
        # Temporary directory
        self.wd = tempfile.mkdtemp()
        # htseq-count fragment
        self.htseq_count = """DDB0166998	2
DDB0167147	0
no_feature	120515
ambiguous	422601
too_low_aQual	0
not_aligned	19996043
alignment_not_unique	0
"""
        # Write to temporary file
        self.htseq_count_file = os.path.join(self.wd,
                                             "htseq_counts.txt")
        with open(self.htseq_count_file,'wt') as fp:
            fp.write(self.htseq_count)
        # Output file name
        self.out_file = os.path.join(self.wd,"out.txt")

    def tearDown(self):
        if os.path.exists(self.wd):
            shutil.rmtree(self.wd)

    def test_annotate_htseq_count_data(self):
        """
        annotate_htseq_count_data: annotates from GFF file
        """
        # Load GFF data
        gff = GFFFile("test.gff",StringIO(gff_data))
        lookup = GFFAnnotationLookup(gff)
        # Do annotation
        annotate_htseq_count_data(lookup,
                                  (self.htseq_count_file,),
                                  self.out_file)
        # Check output
        with open(self.out_file,'rt') as fp:
            self.assertEqual(fp.read(),
                             """exon_parent	feature_type_exon_parent	gene_ID	gene_name	chr	start	end	strand	gene_length	locus	description	htseq_counts.txt
DDB0166998	mRNA	DDB_G0276345	naa20	DDB0232429	6679320	6680012	+	692	DDB0232429:6679320-6680012	Description of gene naa20	2
DDB0167147	mRNA	DDB_G0275629	DDB_G0275629	DDB0232429	5954835	5955486	+	651	DDB0232429:5954835-5955486	Description of gene DDB_G0275629	0
""")
        pass
