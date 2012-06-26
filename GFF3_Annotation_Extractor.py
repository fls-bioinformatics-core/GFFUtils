#!/usr/bin/env python
#
#     GFF3_Annotation_Extractor.py: annotate feature counts with data from GFF
#     Copyright (C) University of Manchester 2012 Peter Briggs, Leo Zeef
#
#######################################################################
#
# GFF3_Annotation_Extractor.py
#
#######################################################################

"""GFF3_Annotation_Extractor.py

Annotate HTSeq-count-style feature count data with information from GFF file.

Usage
-----

First run the htseq-count program to generate counts using e.g.

   htseq-count -q -i ID -t exon <sam> <gff> > htseq-counts

Then run the annotator:

   GFF3_Annotation_Extractor.py -t exon <gff> htseq-counts

Multiple parents
----------------

It's possible for features in a GFF file to have multiple parents.

In this case the output from htseq-count will reproduce the 'Parent'
attribute verbatim, e.g. AF2312,AB2812,abc-3. However
GFF3_Annotation_Extractor will be unable to determine the parent genes
in this case, and so will issue a warning and continue.
"""

#######################################################################
# Module metadata
#######################################################################

import version
__version__ = version.__version__

#######################################################################
# Import modules that this module depends on
#######################################################################

import optparse
import GFFFile
import TabFile
import sys
import logging
import os
import glob

#######################################################################
# Class definitions
#######################################################################

class GFFAnnotationLookup:
    """Utility class for acquiring parent gene names and data

    The GFFAnnotationLookup class provides functionality for indexing
    data from a GFF file in order to facilitate finding the parent
    genes and associated data from the IDs of "feature parents".
    """

    def __init__(self,gff_data):
        """Create a new GFFAnnotationLookup instance

        Arguments:
          gff_data: a GFFFile.GFFFile object populated from a GFF file
        """
        self.__lookup_id = {}
        self.__lookup_parent = {}
        for line in gff_data:
            if 'ID' in line['attributes']:
                # Store reference to data by ID
                idx = line['attributes']['ID']
                self.__lookup_id[idx] = line
                if 'Parent' in line['attributes']:
                    # Store reference to parent by ID
                    parent = line['attributes']['Parent']
                    self.__lookup_parent[idx] = parent

    def getDataFromID(self,idx):
        """Return line of data from GFF file matching the ID attribute

        Arguments:
          idx: ID attribute to search for

        Returns:
          Line of data where the value of the ID attribute matches the
          one supplied; raises KeyError exception if no match is found.
        """
        return self.__lookup_id[idx]

    def getParentID(self,idx):
        """Return ID attribute value for parent feature

        Arguments:
          idx: ID attribute of feature to find the parent of

        Returns:
          ID attribute value of the parent feature; raises KeyError
          exception if no Parent is found
        """
        return self.__lookup_parent[idx]

    def getAncestorGene(self,idx):
        """Return line of data for 'ancestor gene' of feature

        Arguments:
          idx: ID attribute of feature to find the ancestor gene of

        Returns:
          Line of data for the gene which is the ancestor of the
          feature identified by the supplied ID attribute; returns None
          if no parent is found.
        """
        # Follow parents until we find the "ancestor" gene
        idx0 = idx
        try:
            while True:
                idx0 = self.getParentID(idx0)
        except KeyError,ex:
            if idx0 != idx:
                # Check that it's a gene
                data = self.getDataFromID(idx0)
                assert(data['feature'] == 'gene')
                return data
        return None

    def getAnnotation(self,idx):
        """Return annotation data for the supplied feature ID

        Arguments:
          idx: ID attribute of feature to get annotation data for

        Returns:
          GFFAnnotation object populated with the annotation data
          for the feature identified by the supplied ID attribute.
        """
        # Return annotation for an ID
        annotation = GFFAnnotation()
        # Parent feature data
        annotation.parent_feature_name = idx
        parent_feature = self.getDataFromID(idx)
        annotation.parent_feature_type = parent_feature['feature']
        annotation.parent_feature_parent = parent_feature['attributes']['Parent']
        # Parent gene data
        gene = self.getAncestorGene(idx)
        if not gene: return annotation
        annotation.parent_gene_name = gene['attributes']['Name']
        annotation.chr = gene['seqname']
        annotation.start = gene['start']
        annotation.end = gene['end']
        annotation.strand = gene['strand']
        # Locus: chromosome plus start and end data
        annotation.gene_locus = "%s:%s-%s" % (gene['seqname'],gene['start'],gene['end'])
        # Gene length
        annotation.gene_length = gene['end'] - gene['start']
        # Build description text
        # This is all attribute data from the 'description' attribute onwards
        # (but not including the leading "description=" keyword)
        store_attribute = False
        description = []
        for attr in gene['attributes']:
            if store_attribute:
                description.append(attr+'='+gene['attributes'][attr])
            if attr == 'description':
                description.append(gene['attributes'][attr])
                store_attribute = True
        # Reconstruct the description string
        description = ';'.join(description)
        # Finally: replace any tab characters that were introduced by % decoding
        description = description.replace('\t','    ')
        annotation.description = description
        # Done
        return annotation

class GFFAnnotation:
    """Container class for GFF annotation data

    Once instantiated the calling subprogram should populate the
    object properties with the appropriate data.
    """

    def __init__(self):
        """Create a new GFFAnnotation instance
        """
        self.parent_feature_name = ''
        self.parent_feature_type = ''
        self.parent_feature_parent = ''
        self.parent_gene_name = ''
        self.gene_locus = ''
        self.description = ''
        self.chr = ''
        self.start = ''
        self.strand = ''
        self.end = ''
        self.gene_length = ''

#######################################################################
# Functions
#######################################################################

def main():
    """Main program
    """
    # Process command line
    p = optparse.OptionParser(usage="%prog OPTIONS gff_file FEATURE_COUNTS [ FEATURE_COUNTS ... ]",
                              version="%prog "+__version__,
                              description="Annotate feature count data with information from a "
                              "GFF file. Input FEATURE_COUNTS files can be generated using e.g. "
                              "htseq-count ('htseq-count -q -t exon -i Parent gff_file "
                              "sam_file'). The annotator looks up the parent genes of each "
                              "feature and outputs this information against the feature counts "
                              "(in <gff_file>_counts.txt) plus the totals assigned, not "
                              "counted etc (in <gff_file>_counts_stats.txt).")
    p.add_option('-o',action="store",dest="out_file",default=None,
                 help="specify output file name")
    p.add_option('-t','--type',action="store",dest="feature_type",default='exon',
                 help="feature type listed in input count files (default 'exon')")
    options,arguments = p.parse_args()
    if len(arguments) < 2:
        p.error("Expected GFF file and at least one feature count file")

    # Input GFF file
    gff_file = arguments[0]
    if not os.path.exists(gff_file):
        p.error("Input GFF file %s not found" % gff_file)

    # Check for wildcards in feature count file names, to emulate linux shell globbing
    # on platforms such as Windows which don't have this built in
    htseq_files = []
    for arg in arguments[1:]:
        for filen in glob.iglob(arg):
            if not os.path.exists(filen):
                p.error("File '%s' not found" % filen)
            htseq_files.append(filen)
    if not htseq_files:
        p.error("No input feature count files found")

    # Feature type being considered
    feature_type = options.feature_type

    # Output files
    if options.out_file:
        annotated_counts_out_file = options.out_file
    else:
        annotated_counts_out_file = os.path.splitext(os.path.basename(gff_file))[0]+\
            "_counts.txt"
    tables_out_file = \
        os.path.splitext(os.path.basename(annotated_counts_out_file))[0]+\
        "_stats"+os.path.splitext(annotated_counts_out_file)[1]

    # Process GFF data
    print "Reading data from %s" % gff_file
    gff = GFFFile.GFFFile(gff_file)

    # Build lookup
    print "Creating lookup for GFF data"
    gff_lookup = GFFAnnotationLookup(gff)

    # Process the HTSeq-count files
    print "Processing feature count files"
    htseq_counts = GFFFile.OrderedDictionary()
    htseq_tables = GFFFile.OrderedDictionary()
    for htseqfile in htseq_files:
        print "\t%s" % htseqfile
        # Create dictionaries to store data
        htseq_counts[htseqfile] = GFFFile.OrderedDictionary()
        htseq_tables[htseqfile] = GFFFile.OrderedDictionary()
        # Read in data from file
        fp = open(htseqfile,'rU')
        # Flag indicating whether we're reading feature counts
        # or trailing totals
        reading_feature_counts = True
        # Go through the file line-by-line
        for line in fp:
            # All lines are two tab-delimited fields
            name = line.split('\t')[0]
            count = line.strip('\n').split('\t')[1]
            # Check if we've encountered the trailing table
            if line.startswith('no_feature'):
                reading_feature_counts = False
            # Determine what type of data we're storing
            if reading_feature_counts:
                # Feature-by-feature counts i.e.:
                # DDB0166998	1
                htseq_counts[htseqfile][name] = count
            else:
                # Trailing table i.e.:
                # too_low_aQual	0
                htseq_tables[htseqfile][name] = count

    # Total reads counted into genes for each HTSeq output file
    print "Determining numbers counted into genes"
    for htseqfile in htseq_files:
        total_counted_into_genes = 0
        for feature in htseq_counts[htseqfile]:
            total_counted_into_genes += int(htseq_counts[htseqfile][feature])
        # Insert at the start of the table of counts
        htseq_tables[htseqfile].insert(0,
                                       'total_counted_into_genes',
                                       total_counted_into_genes)

    # Make a list of names found in the first count file
    print "Building list of IDs from %s" % htseq_files[0]
    feature_IDs = []
    for feature_ID in htseq_counts[htseq_files[0]]:
        feature_IDs.append(feature_ID)

    # Create a TabFile for output
    print "Building annotated count file for output"
    annotated_counts = TabFile.TabFile(column_names=['exon_parent',
                                                     'feature_type_exon_parent',
                                                     'gene_ID',
                                                     'gene_name',
                                                     'chr',
                                                     'start',
                                                     'end',
                                                     'strand',
                                                     'gene_length',
                                                     'locus',
                                                     'description'])
    for htseqfile in htseq_files:
        annotated_counts.appendColumn(htseqfile)

    # Combine feature counts and parent feature data
    for exon_parent_ID in feature_IDs:
        # Get annotation data
        annotation = gff_lookup.getAnnotation(exon_parent_ID)
        # Build the data line
        data = [annotation.parent_feature_name,
                annotation.parent_feature_type,
                annotation.parent_feature_parent,
                annotation.parent_gene_name,
                annotation.chr,
                annotation.start,
                annotation.end,
                annotation.strand,
                annotation.gene_length,
                annotation.gene_locus,
                annotation.description]
        # Add the counts from each file
        for htseqfile in htseq_files:
            data.append(htseq_counts[htseqfile][exon_parent_ID])
        # Add to the tabfile
        annotated_counts.append(data=data)

    # Write the file
    print "Writing output file %s" % annotated_counts_out_file
    annotated_counts.write(annotated_counts_out_file,include_header=True,no_hash=True)

    # Make second file for the trailing table data
    print "Building trailing tables data file for output"
    table_counts = TabFile.TabFile(column_names=['count'])
    for htseqfile in htseq_files:
        table_counts.appendColumn(htseqfile)
    for name in htseq_tables[htseq_files[0]]:
        # Build the data line
        data = [name]
        for htseqfile in htseq_files:
            data.append(htseq_tables[htseqfile][name])
        table_counts.append(data=data)
    print "Writing output file %s" % tables_out_file
    table_counts.write(tables_out_file,include_header=True,no_hash=True)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s: %(message)s")
    main()
