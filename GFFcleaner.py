#!/bin/env python
#
#     GFFcleaner.py: various functions to clean up GFF files
#     Copyright (C) University of Manchester 2011 Peter Briggs
#
######################################################################
#
# GFFcleaner.py
#
#######################################################################

"""GFFcleaner

Utility program which can perform various operations to 'clean' a GFF
file.

Also includes modules for dealing with GFF files and their 'attribute'
field.
"""
#######################################################################
# Module metadata
#######################################################################

__version__ = "0.2.0"

#######################################################################
# Import modules
#######################################################################

import os,sys
import copy
import logging
import optparse

# Set up for local modules in "share"
SHARE_DIR = os.path.abspath(
    os.path.normpath(
        os.path.join(os.path.dirname(sys.argv[0]),'..','share')))
sys.path.append(SHARE_DIR)
from TabFile import TabFile

#######################################################################
# Classes
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
    name and index properties of the class.

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
# Functions
#######################################################################

def GroupGeneSubsets(gff_data):
    """Group GFF lines into subsets that are the same gene

    Grouping is based on the ID attribute, which is assumed to be of
    the form '<prefix>:<name>:<index>' (e.g. 'CDS:YEL0W:3').

    Consecutive lines are determined to belong to the same gene if
    they have the same name and consecutive indices.

    Arguments:
      gff_data: a list of GFF data lines

    Returns:
      A list of subsets with one or more GFF data lines, where each
      subset is the same gene
    """
    this_subset = []
    subsets = []
    last_index = None
    logging.debug("%d genes submitted for grouping" % len(gff_data))
    for data in gff_data:
        this_index = GFFID(GFFAttributes(data['attributes'])['ID']).index
        if last_index is not None:
            if this_index != last_index + 1:
                subsets.append(this_subset)
                this_subset = []
        this_subset.append(data)
        last_index = this_index
    # Append last subset
    if this_subset:
        subsets.append(this_subset)
    # Report
    for subset in subsets:
        logging.debug("--Subset--")
        for gene in subset:
            logging.debug("\t%s" % GFFID(GFFAttributes(gene['attributes'])['ID']))
    return subsets

def GFFUpdateAttributes(gff_data,update_keys={},exclude_keys=[],no_empty_values=True):
    """Replace and/or exclude data from the GFF attributes

    Performs manipulations on the attribute field of a GFF file, which typically
    consists of key value pairs of the form 'key=value', separated by semicolons.

    The operations are to update the GFF to remove (exclude) specific keys from the
    attributes of all lines, and/or replace the values of keys by copying from
    other keys in the same line.

    Arguments:
      update_keys: a dictionary mapping keys (attribute names) that should be
        replaced with values from other attributes
      exclude_keys: a list of key (attribute names) that should be removed from
        the attribute list
      no_empty_values: if set True (the default) then don't replace existing
        values with blanks (otherwise replacing with blank values is okay)
    """
    for data in gff_data:
        # Process the attributes data
        attributes = GFFAttributes(data['attributes'])
        # Look up values for each of the keys
        for key in attributes.keys():
            # Exclude this key?
            if key in exclude_keys:
                del(attributes[key])
                logging.debug("Excluding %s" % key)
                continue
            # Map to a different value?
            try:
                lookup_key = update_keys[key]
                try:
                    new_value = attributes[lookup_key]
                    if no_empty_values and new_value == '':
                        # If new value is empty then don't replace
                        logging.debug("Not replacing '%s' with empty value '%s'" %
                                      (key,lookup_key))
                    else:
                        attributes[key] = new_value
                except KeyError:
                    logging.warning("No value for '%s' ('%s')" % (lookup_key,key))
            except KeyError:
                # No mapping found for key, ignore
                pass
        # Update attributes in GFF
        data['attributes'] = str(attributes)
        logging.debug("Updated data for output: %s" % data['attributes'])

def GFFGetDuplicateSGDs(gff_data):
    """Return GFF data with duplicate SGD names

    Returns an OrderedDictionary where the keys are the duplicate
    SGD names, e.g.

    dups = GFFGetDuplicateSGDs(gff)
    for sgd in dups.keys():
      ... loops over SGD names ...
      for data in dups[sgd]:
        ... loops over duplicates ...

    Note that duplicates are determined purely by SGD name; no
    account is taken of chromosome or strand.

    Arguments:
      gff_data: a GFFFile object containing the GFF file data

    Returns:
      OrderedDictionary with SGDs as keys for lists of the
      duplicate TabDataLines corresponding to the SGD
    """
    # Use ordered dictionary to sort info on duplicates
    duplicates = OrderedDictionary()
    # Process data line-by-line
    for data in gff_data:
        attributes = GFFAttributes(data['attributes'])
        if 'SGD' in attributes:
            # Store data
            sgd = attributes['SGD']
            if sgd != '':
                # Check for duplicates
                if sgd in duplicates:
                    duplicates[sgd].append(data)
                else:
                    duplicates[sgd] = [data]
    # Filter out true duplicates i.e. SGDs with at least two
    # GFF data lines
    for sgd in duplicates.keys():
            if len(duplicates[sgd]) < 2:
                del(duplicates[sgd])
    # Finished
    return duplicates

def GFFResolveDuplicateSGDs(gff_data,mapping_data,duplicates,overlap_margin):
    """Resolve duplicate SGD names in GFF data

    Attempts to resolve duplicate SGDs by referring to a list of 'best' genes.

    Note that this function doesn't remove any of the data from the GFF input;
    the calling subprogram should do this based on the list of "discards".

    Arguments:

      gff_data: a GFFFile object containing the GFF file data
      mapping_data: a TabFile object containing candidate genes to insert
        into the GFF data if not present
      duplicates: a dictionary with keys representing SGDs (each key maps
        to a list of duplicate GFF data lines for that SGD) returned by the
        GFFGetDuplicateSGDs function
      overlap_margin: additional number of bases either side of candidate
        gene start and end positions to consider when looking for overlaps
        with duplicates

    Returns:

      Dictionary with the following keys (each key linked to a list):

        resolved_sgds: list of SGDs that were resolved
        unresolved_sgds: SGDs that were unresolved
        unresolved_sgds_no_mapping_genes: SGDs that were unresolved
          because there were no mapping genes with the same name
        unresolved_sgds_no_mapping_genes_after_filter: SGDs that were
          unresolved after filtering on chromosome and strand
        unresolved_sgds_no_overlaps: SGDs that were unresolved after
          checking for overlaps
        unresolved_sgds_multiple_matches: SGDs that were unresolved
           because at least two mapping genes matched
        discard: list of GFF data lines which were selected for
           discard after resolution
    """
    # Dictionary for storing results
    result = { 'resolved_sgds': [],
               'unresolved_sgds': [],
               'unresolved_sgds_no_mapping_genes': [],
               'unresolved_sgds_no_mapping_genes_after_filter': [],
               'unresolved_sgds_no_overlaps': [],
               'unresolved_sgds_multiple_matches': [],
               'discard': [] }
    # Make list of matching genes for each duplicate from mapping data
    sgds = duplicates.keys()
    for sgd in duplicates.keys():
        # Look up genes with the same SGD name
        logging.debug("* * * * * * * * * * * * * * * * * * * * * * *")
        logging.debug("SGD = %s" % sgd)
        mapping_genes = mapping_data.lookup('name',sgd)
        if len(mapping_genes) == 0:
            logging.debug("No genes in mapping file with matching SGD to resolve:")
            for duplicate in duplicates[sgd]:
                attr = GFFAttributes(duplicate['attributes'])
                logging.debug("\t%s %s %s %s %s L%d %s" % (attr['ID'],
                                                           duplicate['seqname'],
                                                           duplicate['start'],
                                                           duplicate['end'],
                                                           duplicate['strand'],
                                                           duplicate.lineno(),
                                                           duplicate['feature']))
            # Record unresolved duplicates for reporting
            result['unresolved_sgds_no_mapping_genes'].append(sgd)
            continue
        # At least one mapping gene available
        matches = []
        rejects = []
        # Match duplicates to mapping genes
        genes_to_duplicates = {}
        for duplicate in duplicates[sgd]:
            assigned = False
            for gene in mapping_genes:
                # Filter on chromosome and strand
                if gene['chr'] == duplicate['seqname'] and \
                        gene['strand'] == duplicate['strand']:
                    if gene not in genes_to_duplicates.keys():
                        genes_to_duplicates[gene] = []
                    genes_to_duplicates[gene].append(duplicate)
                    assigned = True
            # No match for this duplicate, add to provisional rejects
            if not assigned:
                if duplicate in rejects:
                    logging.warning("Duplicate added multiple times to rejects list")
                rejects.append(duplicate)
        # Check if there are any matches
        if len(genes_to_duplicates.keys()) == 0:
            logging.debug("No mapping genes matched on chromosome and strand")
            result['unresolved_sgds_no_mapping_genes_after_filter'].append(sgd)
            continue
        # Cluster duplicates for each gene and filter by overlap
        for gene in genes_to_duplicates.keys():
            # Determine overlap region
            region = (gene['start'] - overlap_margin,
                      gene['end'] + overlap_margin)
            # Group duplicates into subsets
            genes_to_duplicates[gene] = GroupGeneSubsets(genes_to_duplicates[gene])
            # Check for overlaps for each subset
            for duplicate in genes_to_duplicates[gene]:
                if duplicate[0]['start'] > region[0] and \
                        duplicate[-1]['end'] < region[1]:
                    # Found a match
                    matches.append(duplicate)
                else:
                    # Not a match, unpack and add to provisional rejects
                    for d in duplicate:
                        if d in rejects:
                            logging.warning("Duplicate added multiple times to rejects list")
                        rejects.append(d)
        # End of filtering process - see what we're left with
        if len(matches) == 1:
            # Resolved
            logging.debug("Duplication resolved for %s" % sgd)
            result['resolved_sgds'].append(sgd)
            for duplicate in matches[0]:
                logging.debug("\t%s %s %s %s L%d %s" % (duplicate['seqname'],
                                                        duplicate['start'],
                                                        duplicate['end'],
                                                        duplicate['strand'],
                                                        duplicate.lineno(),
                                                        duplicate['feature']))
            # Add rejects to discard pile
            for duplicate in rejects:
                result['discard'].append(duplicate)
        elif len(matches) == 0:
            # Unresolved, no overlaps
            result['unresolved_sgds_no_overlaps'].append(sgd)
        else:
            # Multiple matches left
            logging.debug("Unable to resolve duplication for %s between:")
            result['unresolved_sgds_multiple_matches'].append(sgd)
            for match in matches:
                for duplicate in match:
                    logging.debug("\t%s %s %s %s L%d %s" % (duplicate['seqname'],
                                                            duplicate['start'],
                                                            duplicate['end'],
                                                            duplicate['strand'],
                                                            duplicate.lineno(),
                                                            duplicate['feature']))
    # Finished - make list of unresolved SGDs
    for unresolved in ('unresolved_sgds_no_mapping_genes',
                       'unresolved_sgds_no_mapping_genes_after_filter',
                       'unresolved_sgds_no_overlaps',
                       'unresolved_sgds_multiple_matches'):
        result['unresolved_sgds'].extend(result[unresolved])
    return result

def GFFGroupSGDs(gff_data):
    """Update ID attribute of GFF data to indicate SGD groups

    For each line in the GFF data, looks for a non-blank SGD value in the
    GFF attributes and updates the ID attribute with the format:

    ID=CDS:<SGD>:<i>

    where <SGD> is the SGD value and <i> is an integer index starting
    from 1.

    For lines in the GFF within 5 lines of each other and with matching SGDs,
    the integer index increases by 1 each time to indicate that the lines
    form a group, for example CDS:YEL0W:1, CDS:YEL0W:2 etc.

    Arguments:
      gff_data: a GFFFile object containing the GFF file data
    """
    logging.debug("Starting grouping of SGDs")
    next_ln = 0
    for data in gff_data:
        # Increment the line index
        next_ln += 1
        # Process the attributes data
        attributes = GFFAttributes(data['attributes'])
        # Get the SGD value
        sgd = attributes['SGD']
        if sgd != '':
            # Check the ID
            idx = GFFID(attributes['ID'])
            if idx.code != 'CDS':
                # Set the CDS prefix and index and update ID attribute
                idx.code = 'CDS'
                idx.index = 1
                attributes['ID'] = str(idx)
                data['attributes'] = str(attributes)
            ln = data.lineno()
            # Loop over next 5 data lines after this looking for matching SGD
            for data0 in gff_data[next_ln:next_ln+5]:
                attr0 = GFFAttributes(data0['attributes'])
                sgd0 = attr0['SGD']
                if sgd0 == sgd:
                    # Found a match
                    idx0 = GFFID(attr0['ID'])
                    if idx0.code != '':
                        logging.warning("ID already has code assigned (L%d)" % data0.lineno())
                        logging.warning("Index will be overwritten")
                    else:
                        idx0.code = "CDS"
                    idx0.index = idx.index + 1
                    attr0['ID'] = str(idx0)
                    data0['attributes'] = str(attr0)
                    logging.debug("%d %s\t%d %s" % (next_ln,idx,data0.lineno(),idx0))
                    # Don't look any further
                    break
    # Finished grouping by SGD
    return gff_data

def GFFInsertMissingGenes(gff_data,mapping_data):
    """Insert 'missing' genes from mapping file into GFF data

    A gene is considered 'missing' from the GFF data if its name doesn't
    match any of the SGD names in the GFF data.

    Missing genes are inserted into the GFF data at the appropriate
    position based on chromosome and start position.
    
    Arguments:
      gff_data: a GFFFile object containing the GFF file data
      mapping_data: a TabFile object containing candidate genes to insert
        into the GFF data if not present
    """
    # Make a list of all SGDs in current GFF file
    sgds = []
    for data in gff_data:
        attributes = GFFAttributes(data['attributes'])
        if 'SGD' in attributes:
            sgd = attributes['SGD']
            if not sgd in sgds: sgds.append(sgd)
    # Look for SGDs that aren't in the current GFF file
    for gene in mapping_data:
        sgd = gene['name']
        # Throw away anything with invalid data
        try:
            start = int(gene['start'])
        except ValueError:
            logging.warning("Bad start position at L%s ('%s') in %s, skipped" %
                            (gene.lineno(),gene['start'],mapping_data.filename()))
            continue
        if not sgd in sgds:
            # SGD is not in the input GFF
            chrom = gene['chr']
            # Find the insertion point (i.e. within the correct chromosome and
            # start position)
            i = -1
            for j in range(len(gff_data)):
                # Match chromosome
                if gff_data[j]['seqname'] == chrom:
                    if gff_data[j]['start'] < start:
                        i = j + 1
            logging.debug("Inserting '%s' at position %d" % (sgd,i))
            # Insert missing gene into GFF data
            missing = gff_data.insert(i)
            missing['seqname'] = gene['chr']
            missing['source'] = 'GFFcleaner'
            missing['feature'] = 'CDS'
            missing['start'] = gene['start']
            missing['end'] = gene['end']
            missing['score'] = '0'
            missing['strand'] = gene['strand']
            missing['frame'] = '0'
            attributes = GFFAttributes()
            attributes['ID'] = 'CDS:%s:1' % sgd
            attributes['SGD'] = sgd
            attributes['Gene'] = sgd
            attributes['Parent'] = sgd
            missing['attributes'] = str(attributes)
    # Finished inserting missing genes
    return gff_data

#######################################################################
# Tests
#######################################################################

import unittest
import cStringIO

class TestGroupGeneSubsets(unittest.TestCase):

    def setUp(self):
        # List of genes to group
        self.fp = cStringIO.StringIO(
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
        self.fp = cStringIO.StringIO(
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=abc;kaks=-le+100;SGD=YEL0W;ncbi=-1e+100;Name=def;
""")

    def test_update_attributes_exclude_keys(self):
        """Test excluding specific attributes from the attribute field
        """
        gff = GFFFile('test.gff',self.fp)
        # Check that attributes are present initially
        attributes = GFFAttributes(gff[0]['attributes'])
        for attr in ['ID','Name','SGD','ncbi','kaks']:
            self.assertTrue(attr in attributes.keys())
        # Do the exclusion operation
        GFFUpdateAttributes(gff,exclude_keys=['ncbi','kaks'])
        # Check that expected attributes have been removed
        attributes = GFFAttributes(gff[0]['attributes'])
        for attr in ['ID','Name','SGD']:
            self.assertTrue(attr in attributes.keys())
        for attr in ['ncbi','kaks']:
            self.assertTrue(attr not in attributes.keys())

    def test_update_attributes_replace_values(self):
        """Test replacing specific attributes from the attribute field
        """
        gff = GFFFile('test.gff',self.fp)
        GFFUpdateAttributes(gff,update_keys={'ID':'SGD','Name':'SGD'})
        # Check that attributes have the expected values
        attributes = GFFAttributes(gff[0]['attributes'])
        sgd = attributes['SGD']
        for attr in ['ID','Name']:
            self.assertTrue(attributes[attr] == sgd)

class TestGFFGetDuplicateSGDs(unittest.TestCase):

    def setUp(self):
        # Make file-like object to read data in
        self.fp = cStringIO.StringIO(
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
        self.fp = cStringIO.StringIO(
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
        self.mp_resolve_all = cStringIO.StringIO(
"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with a mapping gene removed
        self.mp_missing_mapping_gene = cStringIO.StringIO(
"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t-
""")
        # Mapping data with a mapping gene that doesn't match
        # on strand for one duplicate
        self.mp_missing_matching_mapping_gene = cStringIO.StringIO(
"""YEL0W01\tchr1\t39195\t39569\t-
YEL0W03\tchr1\t34525\t37004\t+
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with a mapping gene that doesn't overlap
        self.mp_mapping_gene_no_overlap = cStringIO.StringIO(
"""YEL0W03\tchr1\t34525\t37004\t-
YEL0W01\tchr1\t41402\t41831\t-
YEL0W02\tchr2\t40406\t40864\t-
""")
        # Mapping data with multiple mapping genes with same name
        self.mp_multiple_mapping_genes = cStringIO.StringIO(
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
        self.fp = cStringIO.StringIO(
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
            idx = GFFAttributes(gff[i]['attributes'])['ID']
            self.assertEqual(idx,self.ids[i],"incorrect ID at position %d" % i)

class TestGFFInsertMissingGenes(unittest.TestCase):

    def setUp(self):
        # Make file-like object for GFF data
        self.fp = cStringIO.StringIO(
"""chr1\tTest\tCDS\t28789\t29049\t0\t-\t0\tID=CDS:YEL0W01:1;SGD=YEL0W01
chr1\tTest\tCDS\t29963\t32155\t0\t-\t0\tID=CDS:YEL0W02:1;SGD=YEL0W02
chr1\tTest\tCDS\t34525\t35262\t0\t-\t0\tID=CDS:YEL0W04:1;SGD=YEL0W04
chr1\tTest\tCDS\t35823\t37004\t0\t-\t0\tID=CDS:YEL0W04:2;SGD=YEL0W04
chr2\tTest\tCDS\t38050\t38120\t0\t-\t0\tID=CDS:YEL0W05:1;SGD=YEL0W05
chr2\tTest\tCDS\t39195\t39569\t0\t-\t0\tID=CDS:YEL0W06:1;SGD=YEL0W06
chr2\tTest\tCDS\t40406\t40864\t0\t-\t0\tID=CDS:YEL0W06:2;SGD=YEL0W06
""")
        # Make a file-like object for mapping data
        self.mp = cStringIO.StringIO(
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
        idx = GFFAttributes(gff[i]['attributes'])['ID']
        self.assertEqual(idx,'CDS:YEL0W03:1',"incorrect ID at position %d" % i)
        # Check: YEL0W06 should not have been inserted
        for i in range(len(gff)):
            idx = GFFAttributes(gff[i]['attributes'])['ID']
            self.assertNotEqual(idx,'CDS:YEL0W06:3',"wrong ID at position %d" %i)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Set up logging format
    logging.basicConfig(format='[%(levelname)s][%(funcName)s] %(message)s')

    p = optparse.OptionParser(usage="%prog [options] <file>.gff",
                              version="%prog "+__version__,
                              description=
                              "Utility to perform various 'cleaning' operations on a GFF file "
                              "and produce output file <file>_clean.gff.")
    p.add_option('--prepend',action='store',dest='prepend_str',default=None,
                 help="String to prepend to seqname in first column")
    p.add_option('--clean',action='store_true',dest='do_clean',
                 help="Perform the 'cleaning' manipulations on the input data")
    p.add_option('--report-duplicates',action='store_true',dest='report_duplicates',
                 help="Report duplicate SGD names and write list to <file>_duplicates.gff "
                 "with line numbers, chromosome, start coordinate and strand.")
    p.add_option('--resolve-duplicates',action='store',dest='mapping_file',default=None,
                 help="Resolve duplicate SGDs by matching against 'best' genes in the supplied "
                 "mapping file; other non-matching genes are discarded and written to "
                 "<file>_discarded.gff.")
    p.add_option('--discard-unresolved',action='store_true',dest='discard_unresolved',
                 help="Also discard any unresolved duplicates, which are written to "
                 "<file>_unresolved.gff.")
    p.add_option('--insert-missing',action='store',dest='gene_file',default=None,
                 help="Insert genes from gene file with SGD names that don't appear in the "
                 "input GFF. If a mapping-file was specified with the --resolve-duplicates "
                 "option then that will be used by default.")
                 
    p.add_option('--debug',action='store_true',dest='debug',
                 help="Print debugging information")
    p.add_option('--test',action='store_true',dest='run_tests',
                 help="Run unit tests")

    # Process the command line
    options,arguments = p.parse_args()

    # Check for debugging
    if options.debug:
        # Turn on debugging output
        logging.getLogger().setLevel(logging.DEBUG)

    # Check for unit testing
    if options.run_tests:
        print "Running unit tests"
        suite = unittest.TestSuite(unittest.TestLoader().\
                                       discover(os.path.dirname(sys.argv[0]), \
                                                    pattern=os.path.basename(sys.argv[0])))
        unittest.TextTestRunner(verbosity=2).run(suite)
        print "Tests finished"
        sys.exit()

    # Input files
    if len(arguments) != 1:
        p.error("input GFF file required")
    else:
        infile = arguments[0]
        if not os.path.exists(infile):
            p.error("Input file '%s' not found" % infile)

    # Report version
    p.print_version()

    # Set flags based on command line

    # String to prepend to first column
    prepend_str = options.prepend_str
    # Cleaning options
    if options.do_clean:
        # Update values in the "score" column
        clean_score = True
        # Clean up the "attributes" column
        clean_attributes = True
        # Initialise mapping of keys from input to output in "attributes" column
        # where new values are required etc
        attributes_key_map = OrderedDictionary()
        attributes_key_map['ID'] = 'SGD'
        attributes_key_map['Gene'] = 'SGD'
        attributes_key_map['Parent'] = 'SGD'
        attributes_key_map['Name'] = 'SGD'
        attributes_dont_replace_with_empty_data = True
        attributes_exclude_keys = ['kaks','kaks2','ncbi']
        # Set ID field in "attributes" to group lines with matching SGDs
        group_SGDs = True
    else:
        clean_score = False
        clean_attributes = False
        group_SGDs = False
    # Report duplicate names
    report_duplicates = options.report_duplicates
    # Resolve duplicated genes using CDS file
    if options.mapping_file is not None:
        resolve_duplicates = True
        cdsfile = options.mapping_file
    else:
        resolve_duplicates = False
        cdsfile = None
    # Discard unresolved duplicates
    discard_unresolved = options.discard_unresolved
    # Insert missing genes
    if options.gene_file is not None:
        insert_missing = True
        if options.gene_file:
            genefile = options.gene_file
        else:
            genefile = cdsfile
    else:
        insert_missing = False
        genefile = None

    # Name for output files
    outext = os.path.splitext(os.path.basename(infile))[1]
    outbase = os.path.splitext(os.path.basename(infile))[0]
    outfile = outbase+'_clean'+outext
    print "Input : %s" % infile
    print "Output: %s" % outfile
    dupfile = outbase+'_duplicates.txt'
    delfile = outbase+'_discarded.gff'
    unresfile = outbase+'_unresolved.gff'

    # Read in data from file
    gff_data = GFFFile(infile)

    # Prepend string to seqname column
    if prepend_str is not None:
        print "Prepending '%s' to values in 'seqname' column" % prepend_str
        for data in gff_data:
            data['seqname'] = prepend_str+str(data['seqname'])

    # Check/clean score column values
    if clean_score:
        print "Replacing 'Anc_*' and blanks with '0's in 'score' column"
        score_unexpected_values = []
        for data in gff_data:
            try:
                # Numerical value
                score = float(data['score'])
                if score != 0:
                    score_unexpected_values.append(data['score'])
            except ValueError:
                # String value
                if data['score'].strip() != '' and not data['score'].startswith('Anc_'):
                    score_unexpected_values.append(data['score'])
            # Replace "Anc_*" or blank values in "score" column with zero
            if data['score'].startswith('Anc_') or data['score'].strip() == '':
                data['score'] = '0'
        # Report unexpected values
        n = len(score_unexpected_values)
        if n > 0:
            logging.warning("%d 'score' values that are not '', 0 or 'Anc_*'" % n)
            logging.warning("Other values: %s" % score_unexpected_values)

    # Clean up the data in "attributes" column
    if clean_attributes:
        print "Cleaning up attributes: replacing keys:"
        for key in attributes_key_map.keys():
            print "\t%s -> %s" % (key,attributes_key_map[key])
        if attributes_dont_replace_with_empty_data:
            print "(Replacement will be skipped if new data is missing/blank)"
        print "Excluding keys:"
        for key in attributes_exclude_keys:
            print "\t%s" % key
        GFFUpdateAttributes(gff_data,attributes_key_map,attributes_exclude_keys,
                            attributes_dont_replace_with_empty_data)

    # Set the IDs for consecutive lines with matching SGD names, to indicate that
    # they're in the same gene
    if group_SGDs:
        print "Grouping SGDs by setting ID's for consecutive lines with the same SGD values"
        GFFGroupSGDs(gff_data)

    # Find duplicates in input file
    if report_duplicates or resolve_duplicates:
        duplicate_sgds = GFFGetDuplicateSGDs(gff_data)
                
    if report_duplicates:
        # Write to duplicates file
        print "Writing duplicate SGD names to %s" % dupfile
        fd = open(dupfile,'w')
        ndup = 0
        ngroups = 0
        for sgd in duplicate_sgds.keys():
            assert(len(duplicate_sgds[sgd]) > 1)
            ndup += 1
            fd.write("%s\t" % sgd)
            for data in duplicate_sgds[sgd]:
                # Write the line number, chromosome, start and strand data
                line = ';'.join(('L'+str(data.lineno()),
                                 str(data['seqname']),str(data['start']),str(data['end'])))
                fd.write("\t%s" % line)
            fd.write("\n")
            logging.debug("%s\t%s" % (sgd,duplicate_sgds[sgd]))
            for group in GroupGeneSubsets(duplicate_sgds[sgd]):
                if len(group) > 1: ngroups += 1
        if ndup == 0:
            fd.write("No duplicate SGDs\n")
        fd.close()
        print "%d duplicates found (of which %d are trivial)" % (ndup,ngroups)

    if resolve_duplicates:
        print "Resolving duplicate SGDs using data from %s" % cdsfile
        print "Discarded genes will be written to %s" % delfile
        # Get data on best gene mappings from CDS file
        # Format is tab-delimited, each line has:
        # orf      chr      start     end      strand
        mapping = TabFile(cdsfile,column_names=('name','chr','start','end','strand'))
        # Overlap margin
        overlap_margin = 1000
        # Perform resolution
        result = GFFResolveDuplicateSGDs(gff_data,mapping,duplicate_sgds,overlap_margin)
        #
        # Report the results
        #
        # Convenience variables for lists of unresolved, discarded etc duplicates
        resolved_sgds = result['resolved_sgds']
        unresolved_sgds_no_mapping_genes = result['unresolved_sgds_no_mapping_genes']
        unresolved_sgds_no_mapping_genes_after_filter = \
            result['unresolved_sgds_no_mapping_genes_after_filter']
        unresolved_sgds_no_overlaps = result['unresolved_sgds_no_overlaps']
        unresolved_sgds_multiple_matches = result['unresolved_sgds_multiple_matches']
        discard = result['discard']
        # Remaining unresolved cases
        if len(unresolved_sgds_no_mapping_genes) > 0:
            print "No mapping genes with same SGDs found in %s:" % cdsfile
            for sgd in unresolved_sgds_no_mapping_genes:
                print "\t%s" % sgd
            print
        if len(unresolved_sgds_no_mapping_genes_after_filter) > 0:
            print "No mapping genes with same chromosome and/or strand:"
            for sgd in unresolved_sgds_no_mapping_genes_after_filter:
                print "\t%s" % sgd
            print
        if len(unresolved_sgds_no_overlaps) > 0:
            print "No mapping genes with overlaps:"
            for sgd in unresolved_sgds_no_overlaps:
                print "\t%s" % sgd
            print
        if len(unresolved_sgds_multiple_matches) > 0:
            print "Multiple matching mapping genes:"
            for sgd in unresolved_sgds_multiple_matches:
                print "\t%s" % sgd
            print
        # Summary counts for each case
        print "Total number of duplicated indexes   : %d" % len(duplicate_sgds.keys())
        print "Number of resolved duplicate SGDs    : %d" % len(resolved_sgds)
        print "Unresolved duplicates:"
        print "* No mapping genes with same SGD     : %d" % len(unresolved_sgds_no_mapping_genes)
        print "* No mapping genes with same chr/str : %d" % len(unresolved_sgds_no_mapping_genes_after_filter)
        print "* No mapping genes with overlap      : %d" % len(unresolved_sgds_no_overlaps)
        print "* Multiple mapping genes match       : %d" % len(unresolved_sgds_multiple_matches)

        # Remove discarded duplicates from the data
        print "Removing discarded duplicates and writing to %s" % delfile
        fd = open(delfile,'w')
        for discard_data in discard:
            try:
                ip = gff_data.indexByLineNumber(discard_data.lineno())
                del(gff_data[ip])
                fd.write("%s\n" % discard_data)
            except IndexError:
                logging.warning("Failed to delete line %d: not found" % discard_data.lineno())
        fd.close()

        # Remove unresolved duplicates if requested
        if discard_unresolved:
            print "Removing unresolved duplicates and writing to %s" % unresfile
            # Get list of unresolved SGDs
            all_unresolved = result['unresolved_sgds']
            # Get list of unresolved duplicates
            unresolved = []
            for data in gff_data:
                attributes = GFFAttributes(data['attributes'])
                if 'SGD' in attributes:
                    if attributes['SGD'] in all_unresolved:
                        unresolved.append(data)
            # Discard them
            fu = open(unresfile,'w')
            for discard in unresolved:
                try:
                    ip = gff_data.indexByLineNumber(discard.lineno())
                    del(gff_data[ip])
                    fu.write("%s\n" % discard)
                except IndexError:
                    logging.warning("Failed to delete line %d: not found" % discard.lineno())
            fu.close()

    # Look for "missing" genes in mapping file
    if insert_missing:
        # Get name for file with gene list
        if genefile is None:
            genefile = cdsfile
        print "Inserting unmatched genes from %s" % genefile
        # Get gene data from CDS file
        # Format is tab-delimited, each line has:
        # orf      chr      start     end      strand
        mapping = TabFile(genefile,column_names=('name','chr','start','end','strand'))
        n_genes_before_insert = len(gff_data)
        gff_data = GFFInsertMissingGenes(gff_data,mapping)
        print "Inserted %d missing genes" % (len(gff_data) - n_genes_before_insert)

    # Write to output file
    print "Writing output file %s" % outfile
    gff_data.write(outfile)


