#!/usr/bin/env python
#
#     annotation.py: handling GFF annotation data
#     Copyright (C) University of Manchester 2020 Peter Briggs
#

import os
import logging
from .GFFFile import OrderedDictionary 
from bcftbx.TabFile import TabFile

#######################################################################
# Class definitions
#######################################################################

class GFFAnnotationLookup(object):
    """Utility class for acquiring parent gene names and data

    The GFFAnnotationLookup class provides functionality for
    indexing data from GFF or GTF files, in order to facilitate
    finding the parent genes and associated data from the IDs of
    "feature parents".

    The following methods are available:

    - getDataFromID: returns the line from the input GFF or GTF
      where the 'ID' attribute (for GFF) or 'gene_id' (for GTF)
      matches the supplied id

    - getParentID: returns the id for the "parent" feature of
      the supplied id, from the 'Parent' attribute (GFF only)

    - getAncestorGene: returns the "ancestor" gene for the
      feature which matches the supplied id; this is determined
      by following the parent features until a 'gene' feature
      is located (GFF only)

    - getAnnotation: returns a GFFAnnotation object for the
      feature that matches the supplied id

    Note that for GTF data only 'gene' features are added to the
    lookup tables.
    """

    def __init__(self,gff_data,id_attr=None,feature_type=None):
        """Create a new GFFAnnotationLookup instance

        Arguments:
          gff_data: a GFFFile object populated from a GFF file
          id_attr: the attribute to use to extract the ID of
            of a feature (defaults to 'ID' for GFF and 'gene_id'
            for GTF, if not set)
          feature_type: if not None then only allow features
            of the specified type to be considered as parents
            when fetching annotation (GFF only)

        """
        self.__feature_data_format = gff_data.format
        self.__lookup_id = {}
        self.__lookup_parent = {}
        self.__feature_type = feature_type
        print("Input file is '%s' format" % self.__feature_data_format)
        if self.__feature_data_format == 'gff':
            self._load_from_gff(gff_data,id_attr=id_attr)
        elif self.__feature_data_format == 'gtf':
            self._load_from_gtf(gff_data,id_attr=id_attr)
        else:
            raise Exception("Unknown format for feature data: '%s'" %
                            gff_data.format)

    def _load_from_gff(self,gff_data,id_attr=None):
        """Create the lookup tables from GFF input
        """
        if id_attr is None:
            id_attr = 'ID'
        parent_attr = 'Parent'
        for line in gff_data:
            if id_attr in line['attributes']:
                # Check that the ID is unique
                idx = line['attributes'][id_attr]
                if idx not in self.__lookup_id:
                    self.__lookup_id[idx] = []
                # Store reference to data by ID
                self.__lookup_id[idx].append(line)
                if parent_attr in line['attributes']:
                    # Store reference to parent by ID
                    parent = line['attributes'][parent_attr]
                    self.__lookup_parent[idx] = parent
                    # Check for multiple parents
                    if len(parent.split(',')) > 1:
                        # Issue a warning but continue for now
                        logging.warning("Multiple parents found on "
                                        "line %d: %s" % (line.lineno(),
                                                         parent))
            else:
                logging.warning("No identifier attribute (%s) on line %d" % 
                                (id_attr,line.lineno()))

    def _load_from_gtf(self,gtf_data,id_attr=None):
        """Create the lookup tables from GTF input
        """
        if id_attr is None:
            id_attr = 'gene_id'
        for line in gtf_data:
            # Only interested in 'gene' features
            if line['feature'] == 'gene':
                if id_attr in line['attributes']:
                    idx = line['attributes'][id_attr]
                    self.__lookup_id[idx] = [line]
                else:
                    logging.warning("No '%s' attribute found on "
                                    "line %d: %s" % (id_attr,
                                                     line.lineno(),
                                                     line))

    def getDataFromID(self,idx):
        """Return line of data from GFF file matching the ID attribute

        Arguments:
          idx: ID attribute to search for

        Returns:
          List of data lines where the value of the ID attribute matches
          the one supplied; raises KeyError exception if no match is
          found.
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
        except KeyError as ex:
            if idx0 != idx:
                # Locate gene record in list
                for data in self.getDataFromID(idx0):
                    if data['feature'] == 'gene':
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
        print("Collecting annotation for %s" % idx)
        # Parent feature data
        parent_feature = None
        try:
            for data in self.getDataFromID(idx):
                if self.__feature_type:
                    # Match the feature type
                    if data['feature'] == self.__feature_type:
                        parent_feature = data
                        break
                else:
                    # Return the first feature
                    parent_feature = data
                    break
        except KeyError:
            # No parent data
            pass
        if not parent_feature:
            logging.warning("Unable to locate parent data for feature "
                            "'%s'" % idx)
            return GFFAnnotation(idx)
        # Parent gene data
        if self.__feature_data_format != 'gtf':
            gene = self.getAncestorGene(idx)
            if not gene:
                return GFFAnnotation(idx,parent_feature)
        else:
            gene = parent_feature
        return GFFAnnotation(idx,parent_feature,gene)

class GFFAnnotation(object):
    """Container class for GFF annotation data

    Once instantiated the calling subprogram should populate the
    object properties with the appropriate data.
    """

    def __init__(self,name,feature=None,associated_gene=None):
        """Create a new GFFAnnotation instance

        Arguments:
          name (str): name of the feature that the
            annotation is associated with
          feature (GFFDataLine): optional feature data
            to populate the annotation from
          associated_gene (GFFDataLine): optional gene
            data to populate the annotation from
        """
        self.parent_feature_name = str(name)
        # Set data from the feature, if supplied
        if feature is not None:
            self.parent_feature_type = feature['feature']
            self.parent_feature_parent = feature['attributes']['Parent']
        else:
            self.parent_feature_type = ''
            self.parent_feature_parent = ''
        # Set data from the associated gene, if supplied
        if associated_gene is not None:
            if associated_gene.format == 'gff':
                self.parent_gene_name = associated_gene['attributes']['Name']
            elif associated_gene.format == 'gtf':
                self.parent_gene_name = associated_gene['attributes']['gene_name']
            self.chr = associated_gene['seqname']
            self.start = associated_gene['start']
            self.end = associated_gene['end']
            self.strand = associated_gene['strand']
            self.description = self.build_description_text(
                associated_gene['attributes'])
        else:
            self.parent_gene_name = ''
            self.description = ''
            self.chr = ''
            self.start = ''
            self.strand = ''
            self.end = ''

    @property
    def gene_locus(self):
        """
        Locus is chromosome plus start and end data
        """
        if self.chr and self.start and self.end:
            return "%s:%s-%s" % (self.chr,self.start,self.end)
        return ''

    @property
    def gene_length(self):
        """
        Gene length is 'end' - 'start'
        """
        if self.start and self.end:
            return int(self.end) - int(self.start)
        return ''

    def build_description_text(self,attributes):
        """
        Build description text from gene attributes data

        This is all attribute data from the 'description'
        attribute onwards
        """
        store_attribute = False
        description = []
        for attr in attributes:
            if store_attribute:
                description.append(attr+'='+attributes[attr])
            if attr == 'description':
                description.append(attributes[attr])
                store_attribute = True
        # Reconstruct the description string
        description = ';'.join(description)
        # Finally: replace any tab characters that were introduced
        # by % decoding
        return description.replace('\t','    ')

class HTSeqCountFile(object):
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
        self.__htseq_counts = OrderedDictionary()
        self.__htseq_table = OrderedDictionary()
        # Total reads counted
        self.__total_reads = 0
        # Read in data from file
        fp = open(htseqfile,'rt')
        # Flag indicating whether we're reading feature counts
        # or trailing totals
        reading_feature_counts = True
        # Go through the file line-by-line
        for line in fp:
            # All lines are two tab-delimited fields
            name = line.split('\t')[0]
            count = line.strip('\n').split('\t')[1]
            # Check if we've encountered the trailing table
            if line.startswith('no_feature') or line.startswith('__no_feature'):
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
                self.__htseq_table[name] = int(count)
        # Finished reading from file
        fp.close()
        # Add total counted at the start of the table of counts
        self.__htseq_table.insert(0,'total_counted_into_genes',
                                  self.__total_reads)

    def feature_IDs(self):
        """Return list of feature IDs from the HTSeq-count output
        """
        return self.__htseq_counts.keys()

    def count(self,feature_id):
        """Return count for feature ID
        """
        return int(self.__htseq_counts[feature_id])

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
    # Determine if input file has a header line
    print("Reading in data from %s" % feature_data_file)
    with open(feature_data_file,'rt') as fp:
        for line in fp:
            if line.startswith('#'):
                first_line_is_header = True
            else:
                first_line_is_header = False
            break

    # Read the feature data into a TabFile
    try:
        input_data = TabFile(filen=feature_data_file,
                             first_line_is_header=first_line_is_header)
    except IndexError as ex:
        if first_line_is_header:
            # Maybe first line was just a comment?
            # Try again without header
            input_data = TabFile(filen=feature_data_file)
        else:
            # Some other failure
            raise ex

    # Initialise columns for output
    if input_data.header():
        # Copy header
        columns = [c for c in input_data.header()]
    else:
        # Make generic header
        columns = ['data%d' % x for x in range(0,input_data.nColumns())]

    # Create and populate second TabFile for output
    feature_data = TabFile(column_names=columns)
    for line in input_data:
        feature_data.append(data=[x for x in line])

    # Append columns for annotation
    print("Appending columns for annotation")
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
    print("Writing output file %s" % out_file)
    feature_data.write(out_file,include_header=True,no_hash=True)

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
        os.path.join(os.path.dirname(annotated_counts_out_file),
                     os.path.splitext(os.path.basename(annotated_counts_out_file))[0]+\
                     "_stats"+os.path.splitext(annotated_counts_out_file)[1])

    # Process the HTSeq-count files
    print("Processing HTSeq-count files")
    htseq_data = {}
    for htseqfile in htseq_files:
        print("\t%s" % htseqfile)
        htseq_data[htseqfile] = HTSeqCountFile(htseqfile)

    # Create a TabFile for output
    print("Building annotated count file for output")
    annotated_counts = TabFile(column_names=['exon_parent',
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
        annotated_counts.appendColumn(os.path.basename(htseqfile))

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
    print("Writing output file %s" % annotated_counts_out_file)
    annotated_counts.write(annotated_counts_out_file,include_header=True,no_hash=True)

    # Make second file for the trailing table data
    print("Building trailing tables data file for output")
    table_counts = TabFile(column_names=['count'])
    for htseqfile in htseq_files:
        table_counts.appendColumn(htseqfile)
    for name in htseq_data[htseq_files[0]].table():
        # Build the data line
        data = [name]
        for htseqfile in htseq_files:
            data.append(htseq_data[htseqfile].table()[name])
        table_counts.append(data=data)
    print("Writing output file %s" % tables_out_file)
    table_counts.write(tables_out_file,include_header=True,no_hash=True)
