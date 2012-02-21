#!/bin/env python
#
#     GFFFile.py: classes for reading and manipulating data in GFF files
#     Copyright (C) University of Manchester 2012 Peter Briggs
#
########################################################################
#
# GFFFile.py
#
#########################################################################

__version__ = "0.1.1"

"""GFFFile

Classes for reading data from GFF (Gene-Finding Format/General Feature
Format) files (see http://www.sanger.ac.uk/resources/software/gff/spec.html
for details of the format).

Classes
-------

There are three GFF-specific classes:

 * GFFFile: read data from GFF so it can be easily interrogated
 * GFFAttributes: read data from GFF attributes field to make it easier to
   handle
 * GFFID: handle data stored in 'ID' attribute

There is an additional base class:

 * OrderedDictionary: augumented dictionary which keeps its keys in the
   order they were added to the dictionary

Usage examples
--------------

Start by reading the data from a GFF file into a GFFFile object:

>>> gff = GFFFile('my.gff')

To iterate over all lines in the GFF and print the feature type:

>>> for line in gff:
>>>    print line['feature']

To iterate over all lines and extract 'ID' field from attributes:

>>> for line in gff:
>>>    attr = GFFAttributes(line['attributes'])
>>>    if 'ID' in attr:
>>>       print attr['ID']

To iterate over all lines and break up 'ID' field into components:

>>> for line in gff:
>>>    attr = GFFAttributes(line['attributes'])
>>>    if 'ID' in attr:
>>>       print attr['ID']

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

from TabFile import TabFile
import logging
import copy

#######################################################################
# Class definitions
#######################################################################

class GFFFile(TabFile):
    """Class for reading GFF files

    See http://www.sanger.ac.uk/resources/software/gff/spec.html
    """
    def __init__(self,gff_file,fp=None,skip_first_line=False,
                 first_line_is_header=False):
        TabFile.__init__(self,gff_file,column_names=('seqname',
                                                     'source',
                                                     'feature',
                                                     'start',
                                                     'end',
                                                     'score',
                                                     'strand',
                                                     'frame',
                                                     'attributes'),
                         fp=fp,skip_first_line=skip_first_line,
                         first_line_is_header=first_line_is_header)

class OrderedDictionary:
    """Augumented dictionary which keeps keys in order

    OrderedDictionary provides an augmented Python dictionary
    class which keeps the dictionary keys in the order they are
    added to the object.

    Items are added, modified and removed as with a standard
    dictionary e.g.:

    >>> d[key] = value
    >>> value = d[key]
    >>> del(d[key])

    The 'keys()' method returns the OrderedDictionary's keys in
    the correct order.
    """
    def __init__(self):
        self.__keys = []
        self.__dict = {}

    def __getitem__(self,key):
        if key not in self.__keys:
            raise KeyError
        return self.__dict[key]

    def __setitem__(self,key,value):
        if key not in self.__keys:
            self.__keys.append(key)
        self.__dict[key] = value

    def __delitem__(self,key):
        try:
            i = self.__keys.index(key)
            del(self.__keys[i])
            del(self.__dict[key])
        except ValueError:
            raise KeyError

    def __len__(self):
        return len(self.__keys)

    def __contains__(self,key):
        return key in self.__keys

    def __iter__(self):
        return iter(self.__keys)

    def keys(self):
        return copy.copy(self.__keys)

class GFFAttributes(OrderedDictionary):
    """Class for handling GFF 'attribute' data

    The GFF 'attribute' data consists of semi-colon separated
    data items, which might be single values or 'key=value'
    pairs.

    Given an input string of this form, a GFFAttributes
    object provides access to the keys via the keys() method
    (lists all keys in the order they were encountered), and
    access to values via the GFFAttributes[key] syntax.

    Keyed data can be modified and deleted via this syntax.

    The values without keys can be accessed as a list via the
    nokeys() method.

    str(GFFAttributes) returns the attribute data as a
    string representation which restores the original format
    found in the GFF file (i.e. semi-column separated data
    items with single values or 'key=value' pairs as
    appropriate).
    """
    def __init__(self,attribute_data=None):
        OrderedDictionary.__init__(self)
        self.__nokeys = []
        # Extract individual data items
        if attribute_data is not None:
            for item in attribute_data.split(';'):
                try:
                    i = item.index('=')
                    key = item[:i].strip()
                    value = item[i+1:].strip()
                except ValueError:
                    # No equals sign
                    key = ''
                    value = item.strip()
                # Store data
                if key == '':
                    # No key: store in a list
                    self.__nokeys.append(value)
                else:
                    # Store key-value pair
                    self[key] = value

    def nokeys(self):
        return self.__nokeys

    def __repr__(self):
        items = []
        for item in self.__nokeys:
            items.append(item)
        for key in self.keys():
            items.append("%s=%s" % (key,self[key]))
        return ';'.join(items)

class GFFID:
    """Class for handling ID attribute in GFF data

    The GFFID class takes an ID string which is expected to be of
    one of the forms

    <name> (e.g. 'XYZ123-A')

    or
    
    <code>:<name>:<index> (e.g. 'CDS:XYZ123-A:2')

    The components can be queried and changed via the code,
    name and index properties of the class, e.g.

    >>> code = GFFID('CDS:XYZ123-A:2').code

    etc

    If there is no code then this is returned as an empty string;
    if there is no index then this is returned as 0.

    Doing e.g. str(gffid) returns the appropriate string
    representation.
    """
    def __init__(self,gff_id):
        items = str(gff_id).split(':')
        if len(items) == 1:
            self.code = ''
            self.name = gff_id
        else:
            self.code = items[0]
            self.name = items[1]
        if len(items) > 2:
            try:
                self.index = int(items[2])
            except ValueError:
                print "Bad id? %s %s" % (gff_id,items)
                raise ValueError
        else:
            self.index = 0
            
    def __repr__(self):
        if not self.code:
            return "%s" % self.name
        else:
            return "%s:%s:%d" % (self.code,self.name,self.index)

#######################################################################
# Tests
#######################################################################

import unittest
import cStringIO

class TestGFFFile(unittest.TestCase):
    """Unit tests for the GFFFile class
    """

    def setUp(self):
        # Example GFF file fragment
        self.fp = cStringIO.StringIO(
"""##gff-version   3
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

class TestOrderDictionary(unittest.TestCase):
    """Unit tests for the OrderedDictionary class
    """

    def test_get_and_set(self):
        """Add and retrieve data
        """
        d = OrderedDictionary()
        self.assertEqual(len(d),0)
        d['hello'] = 'goodbye'
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

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Turn off most logging output for tests
    logging.getLogger().setLevel(logging.CRITICAL)
    # Run tests
    unittest.main()
