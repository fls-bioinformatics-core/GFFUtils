#!/usr/bin/env python

import unittest
from io import StringIO
from bcftbx.TabFile import TabFile
from GFFUtils.GFFFile import GFFFile
from GFFUtils.clean.sgd import *

class TestGroupByID(unittest.TestCase):

    def test_group_by_id(self):
        """
        GroupByID: returns features grouped by ID
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""5	Scannell_and Zill 2011	CDS	3340	5223	Anc_8.272	+	0	ID=CDS:YEL065W:1;SGD=YEL065W;Gene=Sbay_5.2;Parent=Sbay_5.2
5	Scannell_and Zill 2011	CDS	28789	29049	Anc_5.522	-	0	ID=CDS:YDR418W:1;SGD=YDR418W;Gene=Sbay_5.14;Parent=Sbay_5.14
5	Scannell_and Zill 2011	CDS	28789	29313	Anc_5.522	-	0	ID=CDS:YDR418W:2;SGD=YDR418W;Gene=Sbay_5.15;Parent=Sbay_5.15
2	Scannell_and Zill 2011	CDS	1053165	1053662	Anc_5.522	+	0	ID=CDS:YDR418W:3;SGD=YDR418W;Gene=Sbay_2.593;Parent=Sbay_2.593
11	Scannell_and Zill 2011	CDS	623811	624884	Anc_5.716	-	0	ID=CDS:YKR100C:1;SGD=YKR100C;Gene=Sbay_11.338;Parent=Sbay_11.338
11	Scannell_and Zill 2011	CDS	632510	635710	Anc_8.301	+	0	ID=1690892742067829512;SGD=;Gene=Sbay_11.340;Parent=Sbay_11.340
"""))
        groups = GroupByID(gff)
        self.assertEqual(len(groups),4)
        expected_group_ids = ((u"CDS:YEL065W:1",),
                              (u"CDS:YDR418W:1",
                               u"CDS:YDR418W:2",
                               u"CDS:YDR418W:3",),
                              (u"CDS:YKR100C:1",),
                              (u"1690892742067829512",))
        for group,expected_group in zip(groups,expected_group_ids):
            for idx,expected_idx in zip(group,expected_group):
                self.assertEqual(idx['attributes']['ID'],
                                 expected_idx)

class TestGFFGetDuplicateSGDs(unittest.TestCase):

    def test_gff_get_duplicate_sgds(self):
        """
        GFFGetDuplicateSGDs: fetch GFF records with duplicated 'SGD' values
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""5	Scannell_and Zill 2011	CDS	3340	5223	Anc_8.272	+	0	ID=1690892742066571750;ygobhmm=2.2e-19;none=;manual=-1e+100;SGD=YEL065W;Gene=Sbay_5.2;Parent=Sbay_5.2;kaks2=-1e+100;struct=2;ltr=1e+100;hsp=0.0;hcnf=0;aa=0.0;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_8.272;stop=0;introns=0;BLAST=YEL065W;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=1884;hmm=-1e+100;Name=Sbay_5.2;ty=1e+100
5	Scannell_and Zill 2011	CDS	28789	29049	Anc_5.522	-	0	ID=1690892742066567884;ygobhmm=7.9e-55;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_5.14;Parent=Sbay_5.14;kaks2=-1e+100;struct=1;ltr=0.33;hsp=8e-33;hcnf=0;aa=8e-33;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=261;hmm=-1e+100;Name=Sbay_5.14;ty=1e+100
5	Scannell_and Zill 2011	CDS	28789	29313	Anc_5.522	-	0	ID=1690892742066567879;ygobhmm=1.5e-108;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_5.15;Parent=Sbay_5.15;kaks2=-1e+100;struct=2;ltr=1e+100;hsp=4e-65;hcnf=0;aa=4e-65;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=495;hmm=-1e+100;Name=Sbay_5.15;ty=1e+100
2	Scannell_and Zill 2011	CDS	1053165	1053662	Anc_5.522	+	0	ID=1690892742067239098;ygobhmm=3.2e-111;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_2.593;Parent=Sbay_2.593;kaks2=-1e+100;struct=2;ltr=0.64;hsp=1e-67;hcnf=0;aa=1e-67;lda=-1e+100;rna=-1e+100;synt=3;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=498;hmm=-1e+100;Name=Sbay_2.593;ty=1e+100
11	Scannell_and Zill 2011	CDS	623811	624884	Anc_5.716	-	0	ID=1690892742067852754;ygobhmm=2.4e-104;none=;manual=-1e+100;SGD=YKR100C;Gene=Sbay_11.338;Parent=Sbay_11.338;kaks2=-1e+100;struct=2;ltr=1.4;hsp=6e-115;hcnf=0;aa=6e-115;lda=-1e+100;rna=-1e+100;synt=4;YGOB=Anc_5.716;stop=0;introns=0;BLAST=YKR100C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=1074;hmm=-1e+100;Name=Sbay_11.338;ty=5.7
"""))
        duplicates = GFFGetDuplicateSGDs(gff)
        self.assertEqual(list(duplicates.keys()),
                         ['YDR418W'])
        expected_duplicates = (
            "5	Scannell_and Zill 2011	CDS	28789	29049	Anc_5.522	-	0	ID=1690892742066567884;ygobhmm=7.9e-55;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_5.14;Parent=Sbay_5.14;kaks2=-1e+100;struct=1;ltr=0.33;hsp=8e-33;hcnf=0;aa=8e-33;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=261;hmm=-1e+100;Name=Sbay_5.14;ty=1e+100",
            "5	Scannell_and Zill 2011	CDS	28789	29313	Anc_5.522	-	0	ID=1690892742066567879;ygobhmm=1.5e-108;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_5.15;Parent=Sbay_5.15;kaks2=-1e+100;struct=2;ltr=1e+100;hsp=4e-65;hcnf=0;aa=4e-65;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=495;hmm=-1e+100;Name=Sbay_5.15;ty=1e+100",
            "2	Scannell_and Zill 2011	CDS	1053165	1053662	Anc_5.522	+	0	ID=1690892742067239098;ygobhmm=3.2e-111;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_2.593;Parent=Sbay_2.593;kaks2=-1e+100;struct=2;ltr=0.64;hsp=1e-67;hcnf=0;aa=1e-67;lda=-1e+100;rna=-1e+100;synt=3;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=498;hmm=-1e+100;Name=Sbay_2.593;ty=1e+100",)
        for dup,expected in zip(duplicates['YDR418W'],
                                expected_duplicates):
            self.assertEqual(str(dup),expected)

class TestGFFResolveDuplicateSGDs(unittest.TestCase):

    def test_gff_resolve_duplicate_sgds(self):
        """
        GFFResolveDuplicateSGDs: resolve duplicate SGDs using mapping file
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""chr5	Scannell_and Zill 2011	CDS	3340	5223	Anc_8.272	+	0	ID=1690892742066571750;ygobhmm=2.2e-19;none=;manual=-1e+100;SGD=YEL065W;Gene=Sbay_5.2;Parent=Sbay_5.2;kaks2=-1e+100;struct=2;ltr=1e+100;hsp=0.0;hcnf=0;aa=0.0;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_8.272;stop=0;introns=0;BLAST=YEL065W;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=1884;hmm=-1e+100;Name=Sbay_5.2;ty=1e+100
chr5	Scannell_and Zill 2011	CDS	28789	29049	Anc_5.522	-	0	ID=1690892742066567884;ygobhmm=7.9e-55;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_5.14;Parent=Sbay_5.14;kaks2=-1e+100;struct=1;ltr=0.33;hsp=8e-33;hcnf=0;aa=8e-33;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=261;hmm=-1e+100;Name=Sbay_5.14;ty=1e+100
chr2	Scannell_and Zill 2011	CDS	1053165	1053662	Anc_5.522	+	0	ID=1690892742067239098;ygobhmm=3.2e-111;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_2.593;Parent=Sbay_2.593;kaks2=-1e+100;struct=2;ltr=0.64;hsp=1e-67;hcnf=0;aa=1e-67;lda=-1e+100;rna=-1e+100;synt=3;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=498;hmm=-1e+100;Name=Sbay_2.593;ty=1e+100
chr11	Scannell_and Zill 2011	CDS	623811	624884	Anc_5.716	-	0	ID=1690892742067852754;ygobhmm=2.4e-104;none=;manual=-1e+100;SGD=YKR100C;Gene=Sbay_11.338;Parent=Sbay_11.338;kaks2=-1e+100;struct=2;ltr=1.4;hsp=6e-115;hcnf=0;aa=6e-115;lda=-1e+100;rna=-1e+100;synt=4;YGOB=Anc_5.716;stop=0;introns=0;BLAST=YKR100C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=1074;hmm=-1e+100;Name=Sbay_11.338;ty=5.7
"""))
        mapping_data = TabFile(fp=StringIO(u"""YEL065W	chr5	271736	353971	+
YDR418W	chr5	28894	29043	-
YKR100C	chr5	1050499	1112150	-
"""),
                               column_names=('name',
                                             'chr',
                                             'start',
                                             'end',
                                             'strand'))
        duplicates = GFFGetDuplicateSGDs(gff)
        results = GFFResolveDuplicateSGDs(gff,
                                          mapping_data,
                                          duplicates)
        self.assertEqual(results['resolved_sgds'],['YDR418W'])
        self.assertEqual(results['unresolved_sgds'],[])
        self.assertEqual(results['unresolved_sgds_no_mapping_genes'],[])
        self.assertEqual(results['unresolved_sgds_no_mapping_genes_after_filter'],[])
        self.assertEqual(results['unresolved_sgds_no_overlaps'],[])
        self.assertEqual(results['unresolved_sgds_multiple_matches'],[])
        self.assertEqual(len(results['discard']),1)

class TestGFFGroupSGDs(unittest.TestCase):

    def test_gff_group_sgds(self):
        """
        GFFGroupSGDs: relabel ID with SGD information
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""5	Scannell_and Zill 2011	CDS	3340	5223	Anc_8.272	+	0	ID=1690892742066571750;ygobhmm=2.2e-19;none=;manual=-1e+100;SGD=YEL065W;Gene=Sbay_5.2;Parent=Sbay_5.2;kaks2=-1e+100;struct=2;ltr=1e+100;hsp=0.0;hcnf=0;aa=0.0;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_8.272;stop=0;introns=0;BLAST=YEL065W;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=1884;hmm=-1e+100;Name=Sbay_5.2;ty=1e+100
5	Scannell_and Zill 2011	CDS	28789	29049	Anc_5.522	-	0	ID=1690892742066567884;ygobhmm=7.9e-55;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_5.14;Parent=Sbay_5.14;kaks2=-1e+100;struct=1;ltr=0.33;hsp=8e-33;hcnf=0;aa=8e-33;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=261;hmm=-1e+100;Name=Sbay_5.14;ty=1e+100
5	Scannell_and Zill 2011	CDS	28789	29313	Anc_5.522	-	0	ID=1690892742066567879;ygobhmm=1.5e-108;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_5.15;Parent=Sbay_5.15;kaks2=-1e+100;struct=2;ltr=1e+100;hsp=4e-65;hcnf=0;aa=4e-65;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=495;hmm=-1e+100;Name=Sbay_5.15;ty=1e+100
2	Scannell_and Zill 2011	CDS	1053165	1053662	Anc_5.522	+	0	ID=1690892742067239098;ygobhmm=3.2e-111;none=;manual=-1e+100;SGD=YDR418W;Gene=Sbay_2.593;Parent=Sbay_2.593;kaks2=-1e+100;struct=2;ltr=0.64;hsp=1e-67;hcnf=0;aa=1e-67;lda=-1e+100;rna=-1e+100;synt=3;YGOB=Anc_5.522;stop=0;introns=0;BLAST=YEL054C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=498;hmm=-1e+100;Name=Sbay_2.593;ty=1e+100
11	Scannell_and Zill 2011	CDS	623811	624884	Anc_5.716	-	0	ID=1690892742067852754;ygobhmm=2.4e-104;none=;manual=-1e+100;SGD=YKR100C;Gene=Sbay_11.338;Parent=Sbay_11.338;kaks2=-1e+100;struct=2;ltr=1.4;hsp=6e-115;hcnf=0;aa=6e-115;lda=-1e+100;rna=-1e+100;synt=4;YGOB=Anc_5.716;stop=0;introns=0;BLAST=YKR100C;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=1074;hmm=-1e+100;Name=Sbay_11.338;ty=5.7
11	Scannell_and Zill 2011	CDS	632510	635710	Anc_8.301	+	0	ID=1690892742067829512;ygobhmm=9e-33;none=;manual=-1e+100;SGD=;Gene=Sbay_11.340;Parent=Sbay_11.340;kaks2=-1e+100;struct=-1;ltr=1.1;hsp=3e-86;hcnf=0;aa=3e-86;lda=-1e+100;rna=-1e+100;synt=0;YGOB=Anc_8.301;stop=0;introns=0;BLAST=YKR102W;ncbi=1e+100;kaks=1e+100;nnnn=-1e+100;length=3201;hmm=-1e+100;Name=Sbay_11.340;ty=1e+100
"""))
        GFFGroupSGDs(gff)
        expected_ids = ("CDS:YEL065W:1",
                        "CDS:YDR418W:1",
                        "CDS:YDR418W:2",
                        "CDS:YDR418W:3",
                        "CDS:YKR100C:1",
                        "1690892742067829512")
        for line,expected_id in zip(gff,expected_ids):
            self.assertEqual(line['attributes']['ID'],expected_id)

class TestGFFInsertMissingGenes(unittest.TestCase):

    def test_gff_insert_missing_genes(self):
        """
        GFFInsertMissingGenes: insert gene from mapping data file
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""chr5	Scannell_and Zill 2011	CDS	3340	5223	Anc_8.272	+	0	ID=1690892742066571750;SGD=YEL065W;Gene=Sbay_5.2;Parent=Sbay_5.2;
chr5	Scannell_and Zill 2011	CDS	28789	29049	Anc_5.522	-	0	ID=1690892742066567884;SGD=YDR418W;Gene=Sbay_5.14;Parent=Sbay_5.14
chr5	Scannell_and Zill 2011	CDS	28789	29313	Anc_5.522	-	0	ID=1690892742066567879;SGD=YDR418W;Gene=Sbay_5.15;Parent=Sbay_5.15
chr2	Scannell_and Zill 2011	CDS	1053165	1053662	Anc_5.522	+	0	ID=1690892742067239098;SGD=YDR418W;Gene=Sbay_2.593;Parent=Sbay_2.593
chr11	Scannell_and Zill 2011	CDS	623811	624884	Anc_5.716	-	0	ID=1690892742067852754;SGD=YKR100C;Gene=Sbay_11.338;Parent=Sbay_11.338
chr11	Scannell_and Zill 2011	CDS	632510	635710	Anc_8.301	+	0	ID=1690892742067829512;SGD=;Gene=Sbay_11.340;Parent=Sbay_11.340
"""))
        mapping_data = TabFile(fp=StringIO(u"""YEL065W	chr5	271736	353971	+
YDR418W	chr5	28894	29043	-
YKR100C	chr5	1050499	1112150	-
YIL151C	chr5	1887446	1900607	+
"""),
                               column_names=('name',
                                             'chr',
                                             'start',
                                             'end',
                                             'strand'))
        GFFInsertMissingGenes(gff,mapping_data)
        expected_gff = GFFFile("expected.gff",
                               fp=StringIO(u"""chr5	Scannell_and Zill 2011	CDS	3340	5223	Anc_8.272	+	0	ID=1690892742066571750;SGD=YEL065W;Gene=Sbay_5.2;Parent=Sbay_5.2;
chr5	Scannell_and Zill 2011	CDS	28789	29049	Anc_5.522	-	0	ID=1690892742066567884;SGD=YDR418W;Gene=Sbay_5.14;Parent=Sbay_5.14
chr5	Scannell_and Zill 2011	CDS	28789	29313	Anc_5.522	-	0	ID=1690892742066567879;SGD=YDR418W;Gene=Sbay_5.15;Parent=Sbay_5.15
chr5	GFFcleaner	CDS	1887446	1900607	0	+	0	ID=CDS:YIL151C:1;SGD=YIL151C;Gene=YIL151C;Parent=YIL151C
chr2	Scannell_and Zill 2011	CDS	1053165	1053662	Anc_5.522	+	0	ID=1690892742067239098;SGD=YDR418W;Gene=Sbay_2.593;Parent=Sbay_2.593
chr11	Scannell_and Zill 2011	CDS	623811	624884	Anc_5.716	-	0	ID=1690892742067852754;SGD=YKR100C;Gene=Sbay_11.338;Parent=Sbay_11.338
chr11	Scannell_and Zill 2011	CDS	632510	635710	Anc_8.301	+	0	ID=1690892742067829512;SGD=;Gene=Sbay_11.340;Parent=Sbay_11.340
"""))
        for line,expected in zip(gff,expected_gff):
            self.assertEqual(str(line),str(expected))
