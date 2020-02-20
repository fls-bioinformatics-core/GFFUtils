#!/usr/bin/env python

import unittest
from io import StringIO
from GFFUtils.GFFFile import GFFFile
from GFFUtils.clean.generic import *

class TestGFFUpdateAttributes(unittest.TestCase):

    def test_gff_update_attributes_defaults(self):
        """
        GFFUpdateAttributes: defaults don't alter lines
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""89	Scannell_and Zill 2011	gap	283	504		+	-	ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100"""))
        GFFUpdateAttributes(gff)
        self.assertEqual(str(gff[0]['attributes']),
                         "ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100")

    def test_gff_update_attributes_update_keys(self):
        """
        GFFUpdateAttributes: update keys with new values
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""89	Scannell_and Zill 2011	gap	283	504		+	-	ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100"""))
        GFFUpdateAttributes(gff,
                            update_keys={
                                'SGD': 'ID',
                                'ty': 'struct'
                            })
        self.assertEqual(str(gff[0]['attributes']),
                         "ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=1690892742066564085;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=0")

    def test_gff_update_attributes_update_allow_empty_value(self):
        """
        GFFUpdateAttributes: update and allow empty value
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""89	Scannell_and Zill 2011	gap	283	504		+	-	ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100"""))
        GFFUpdateAttributes(gff,
                            update_keys={ 'Gene':'SGD' },
                            no_empty_values=False)
        self.assertEqual(str(gff[0]['attributes']),
                         "ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100")

    def test_gff_update_attributes_update_dont_allow_empty_value(self):
        """
        GFFUpdateAttributes: update but don't allow empty value
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""89	Scannell_and Zill 2011	gap	283	504		+	-	ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100"""))
        GFFUpdateAttributes(gff,
                            update_keys={ 'Gene':'SGD' },
                            no_empty_values=True,
        )
        self.assertEqual(str(gff[0]['attributes']),
                         "ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100")

    def test_gff_update_attributes_exclude_keys(self):
        """
        GFFUpdateAttributes: exclude specified keys
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""89	Scannell_and Zill 2011	gap	283	504		+	-	ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100"""))
        GFFUpdateAttributes(gff,exclude_keys=('kaks','kaks2'))
        self.assertEqual(str(gff[0]['attributes']),
                         "ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100")

    def test_gff_update_attributes_exclude_nokeys(self):
        """
        GFFUpdateAttributes: exclude values which don't have associated keys
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""89	Scannell_and Zill 2011	gap	283	504		+	-	ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100;testing;Sbay_89.1"""))
        GFFUpdateAttributes(gff,exclude_nokeys=True)
        self.assertEqual(str(gff[0]['attributes']),
                         "ID=1690892742066564085;ygobhmm=1e+100;none=;manual=-1e+100;SGD=;Gene=Sbay_89.1;Parent=Sbay_89.1;kaks2=-1e+100;struct=0;ltr=1e+100;hsp=1e+100;hcnf=0;aa=1e+100;lda=-1e+100;rna=-1e+100;synt=0;YGOB=;stop=0;introns=0;BLAST=;ncbi=1e+100;kaks=1e+100;nnnn=1e+100;length=222;hmm=-1e+100;Name=Sbay_89.1;ty=1e+100")

class TestGFFAddExonIDs(unittest.TestCase):

    def test_gff_add_exon_ids(self):
        """
        GFFAddExonIDs: constructs and adds missing ID attributes for exons
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""##gff-version   3
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
DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=Description of gene naa20"""))
        GFFAddExonIDs(gff)
        expected_ids = ("DDB0166998",
                        "exon:DDB0166998:00000001",
                        None,
                        "exon:DDB0166998:00000002",
                        None,
                        "exon:DDB0166998:00000003",
                        None,
                        "DDB_G0276345")
        for line,expected_id in zip(gff,expected_ids):
            if expected_id is not None:
                self.assertEqual(line['attributes']['ID'],expected_id)
            else:
                try:
                    line['attributes']['ID']
                    self.fail()
                except KeyError:
                    pass
        

class TestGFFAddIDAttributes(unittest.TestCase):

    def test_gff_add_id_attributes(self):
        """
        GFFAddIDAttributes: constructs and adds missing ID attributes
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""##gff-version   3
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
DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=Description of gene naa20"""))
        GFFAddIDAttributes(gff)
        expected_ids = ("DDB0166998",
                        "exon:DDB0166998:00000001",
                        "CDS:DDB0166998:00000002",
                        "exon:DDB0166998:00000003",
                        "CDS:DDB0166998:00000004",
                        "exon:DDB0166998:00000005",
                        "CDS:DDB0166998:00000006",
                        "DDB_G0276345")
        for line,expected_id in zip(gff,expected_ids):
            self.assertEqual(line['attributes']['ID'],expected_id)

class TestGFFDecodeAttributes(unittest.TestCase):

    def test_gff_decode_attributes(self):
        """
        GFFDecodeAttributes: re-encode escape sequences
        """
        gff = GFFFile("test.gff",
                      fp=StringIO(u"""##gff-version   3
##sequence-region DDB0232428 1 4923596
##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=352472
# generated: Wed Jan 25 03:46:58 2012
DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=ortholog of the conserved catalytic subunit of the NatB N-terminal acetyltransferase (NAA20)%2C which in yeast catalyzes the transfer of an acetyl group to the N-terminal residue of a protein that contains a Met-Glu%2C Met-Asp%2C Met-Asn%2C or Met-Met N-terminus"""))
        GFFDecodeAttributes(gff)
        self.assertEqual(str(gff[0]),
                         "DDB0232429	.	gene	6679320	6680012	.	+	.	ID=DDB_G0276345;Name=naa20;description=ortholog of the conserved catalytic subunit of the NatB N-terminal acetyltransferase (NAA20), which in yeast catalyzes the transfer of an acetyl group to the N-terminal residue of a protein that contains a Met-Glu, Met-Asp, Met-Asn, or Met-Met N-terminus")

