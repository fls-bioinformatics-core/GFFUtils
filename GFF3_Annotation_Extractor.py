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

Annotate gene feature data (for example the output from one or more runs of the
HTSeq-count program) by combining it with data about each feature's parent gene,
taken from a GFF file.

By default the program takes a single tab-delimited input file where the first column
contains feature IDs, and appends data about the feature's parent gene.

In "htseq-count" mode, one or more `htseq-count` output files should be provided as
input, and the program will write out the data about the feature's parent gene appended
with the counts from each input file.

To generate the feature count files using `htseq-count` do e.g.:

        htseq-count --type=exon -i Parent <file>.gff <file>.sam

which returns counts of each exon against the name of that exon's parent.

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
                    # Check for multiple parents
                    if len(parent.split(',')) > 1:
                        # Issue a warning but continue for now
                        logging.warning("Multiple parents found on line %d: %s" % (line.lineno(),
                                                                                   parent))

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
        try:
            parent_feature = self.getDataFromID(idx)
        except KeyError:
            # No parent data
            logging.warning("No parent data for feature '%s'" % idx)
            return annotation
        annotation.parent_feature_name = idx
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

class HTSeqCountFile:
    """Class for handling data from output of htseq-count program

    The htseq-count program outputs 2-columns of tab-delimited data,
    with the first column containing feature IDs and the second the
    corresponding count of reads assigned to that feature. Appended
    to this are summary statistics (for example count of reads not
    assigned to any feature).

    The HTSeqCountFile class processes the data from this output file,
    splitting the data into counts against features and a 'trailing
    table' for the summary statistic data.

    A list of feature IDs can be obtained using the feature_IDs()
    method; the count for a specified feature ID can be obtained
    via the count() method; and the trailing table can be obtained
    using the table() method.
    """

    def __init__(self,htseqfile):
        """Create new HTSeqCountFile instance

        Arguments:
          htseqfile: name of the htseq-count output file (including
            leading path) to process
        """
        # Create dictionaries to store data
        self.__htseq_counts = GFFFile.OrderedDictionary()
        self.__htseq_table = GFFFile.OrderedDictionary()
        # Total reads counted
        self.__total_reads = 0
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
                self.__htseq_counts[name] = count
                self.__total_reads += int(count)
            else:
                # Trailing table i.e.:
                # too_low_aQual	0
                self.__htseq_table[name] = count
        # Finished reading from file
        fp.close()
        # Add total counted at the start of the table of counts
        self.__htseq_table.insert(0,'total_counted_into_genes',
                                  self.__total_reads)

    def feature_IDs(self):
        """Return list of feature IDs from the HTSeq-count output
        """
        return self.__htseq_counts

    def count(self,feature_id):
        """Return count for feature ID
        """
        return self.__htseq_counts[feature_id]

    def table(self):
        """Return the trailing table data

        Returns:
          OrderDictionary object with the statistics as keys referencing
          the values as dictionary items.
        """
        return self.__htseq_table

#######################################################################
# Functions
#######################################################################

# annotate_feature_data
#
def annotate_feature_data(gff_lookup,feature_data_file,out_file):
    """Annotate feature data with gene information

    Reads in 'feature data' from a tab-delimited input file with feature
    IDs in the first column; outputs these data with data about the
    parent gene appended to each line.

    Arguments:
      gff_lookup         populated GFFAnnotationLookup instance
      feature_data_file  input data file with feature IDs in first column
      out_file           name of output file
    """
    # Read the feature data into a TabFile
    print "Reading in data from %s" % feature_data_file
    feature_data = TabFile.TabFile(filen=feature_data_file,
                                   first_line_is_header=True)

    # Append columns for annotation
    print "Appending columns for annotation"
    for colname in ('exon_parent',
                    'feature_type_exon_parent',
                    'gene_ID',
                    'gene_name',
                    'chr',
                    'start',
                    'end',
                    'strand',
                    'gene_length',
                    'locus',
                    'description'):
        feature_data.appendColumn(colname)

    for line in feature_data:
        feature_ID = line[0]
        annotation = gff_lookup.getAnnotation(feature_ID)
        line['exon_parent'] = annotation.parent_feature_name
        line['feature_type_exon_parent'] = annotation.parent_feature_type
        line['gene_ID'] = annotation.parent_feature_parent
        line['gene_name'] = annotation.parent_gene_name
        line['chr'] = annotation.chr
        line['start'] = annotation.start
        line['end'] = annotation.end
        line['strand'] = annotation.strand
        line['gene_length'] = annotation.gene_length
        line['locus'] = annotation.gene_locus
        line['description'] = annotation.description

    # Output
    print "Writing output file %s" % out_file
    feature_data.write(out_file,include_header=True,no_hash=True)

# annotate_htseq_count_data
#
def annotate_htseq_count_data(gff_lookup,htseq_files,out_file):
    """Annotate count data from htseq-count output with gene information

    Reads in data from one or more htseq-count output files and combines
    into a single tab-delimited output file where the counts for each
    feature have been appended to data about the parent gene.

    Also creates an output 'stats' file which combines the summary data
    from the tail of each htseq-count file.

    Arguments:
      gff_lookup:  populated GFFAnnotationLookup instance
      htseq_files: list of output files from htseq-count to use as input
      out_file:    name of output file
    """
    # Output files
    annotated_counts_out_file = out_file
    tables_out_file = \
        os.path.splitext(os.path.basename(annotated_counts_out_file))[0]+\
        "_stats"+os.path.splitext(annotated_counts_out_file)[1]

    # Process the HTSeq-count files
    print "Processing HTSeq-count files"
    htseq_data = {}
    for htseqfile in htseq_files:
        print "\t%s" % htseqfile
        htseq_data[htseqfile] = HTSeqCountFile(htseqfile)

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
    for feature_ID in htseq_data[htseq_files[0]].feature_IDs():
        # Get annotation data
        annotation = gff_lookup.getAnnotation(feature_ID)
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
            data.append(htseq_data[htseqfile].count(feature_ID))
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
    for name in htseq_data[htseq_files[0]].table():
        # Build the data line
        data = [name]
        for htseqfile in htseq_files:
            data.append(htseq_data[htseqfile].table()[name])
        table_counts.append(data=data)
    print "Writing output file %s" % tables_out_file
    table_counts.write(tables_out_file,include_header=True,no_hash=True)

# Main program
#
def main():
    """Main program
    """
    # Process command line
    p = optparse.OptionParser(usage="\n  %prog OPTIONS GFF_FILE FEATURE_DATA\n"
                              "  %prog --htseq-count OPTIONS GFF_FILE FEATURE_COUNTS "
                              "[ FEATURE_COUNTS ... ]",
                              version="%prog "+__version__,
                              description="Annotate feature count data with information from a "
                              "GFF file. Default mode is to take a single tab-delimited "
                              "FEATURE_DATA input file where the first column consists of feature "
                              "IDs from the input GFF_FILE; in this mode each line of "
                              "FEATURE_DATA will be appended with data about the 'parent feature' "
                              "and 'parent gene' matching the feature ID. In --htseq-count mode "
                              "input consists of one or more FEATURE_COUNTS files generated using "
                              "htseq-count (e.g. 'htseq-count -q -t exon -i Parent gff_file "
                              "sam_file'). The annotator looks up the parent genes of each "
                              "feature and outputs this information against the feature counts "
                              "(in <gff_file>_annot.txt) plus the totals assigned, not "
                              "counted etc (in <gff_file>_annot_stats.txt).")
    p.add_option('-o',action="store",dest="out_file",default=None,
                 help="specify output file name")
    p.add_option('-t','--type',action="store",dest="feature_type",default='exon',
                 help="feature type listed in input count files (default 'exon')")
    p.add_option('--htseq-count',action="store_true",dest="htseq_count",default=False,
                 help="htseq-count mode: input is one or more output FEATURE_COUNT files from "
                 "the htseq-count program")
    options,arguments = p.parse_args()

    # Determine what mode to operate in
    htseq_count_mode = options.htseq_count

    # Initial check on arguments
    if len(arguments) < 2:
        p.error("Expected GFF file and at least one feature data file")

    # Input GFF file
    gff_file = arguments[0]
    if not os.path.exists(gff_file):
        p.error("Input GFF file %s not found" % gff_file)

    # Check for wildcards in feature data file names, to emulate linux shell globbing
    # on platforms such as Windows which don't have this built in
    feature_data_files = []
    for arg in arguments[1:]:
        for filen in glob.iglob(arg):
            if not os.path.exists(filen):
                p.error("File '%s' not found" % filen)
            feature_data_files.append(filen)
    if not feature_data_files:
        p.error("No input feature data files found")

    # Final check on number of input files
    if not htseq_count_mode and len(feature_data_files) > 1:  
        p.error("Expected GFF file and a single feature data file")

    # Feature type being considered
    feature_type = options.feature_type

    # Output file
    if options.out_file:
        out_file = options.out_file
    else:
        out_file = os.path.splitext(os.path.basename(gff_file))[0] + "_annot.txt"

    # Process GFF data
    print "Reading data from %s" % gff_file
    gff = GFFFile.GFFFile(gff_file)

    # Build lookup
    print "Creating lookup for GFF data"
    gff_lookup = GFFAnnotationLookup(gff)

    # Annotate input data
    if htseq_count_mode:
        # HTSeq-count mode
        annotate_htseq_count_data(gff_lookup,
                                  feature_data_files,
                                  out_file)
    else:
        # Standard mode
        annotate_feature_data(gff_lookup,
                              feature_data_files[0],
                              out_file)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s: %(message)s")
    main()
