#!/usr/bin/env python
#
#     GFF_HTSeq_Annotator.py: annotate HTSeq-count output with data from GFF
#     Copyright (C) University of Manchester 2012 Peter Briggs, Leo Zeef
#
#######################################################################
#
# GFF_HTSeq_Annotator.py
#
#######################################################################

"""GFF_HTSeq_Annotator.py

Annotate HTSeq-count output with data from GFF
"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.0.1"

#######################################################################
# Import modules that this module depends on
#######################################################################

import optparse
import GFFcleaner
import TabFile
import sys
import logging
import os

#######################################################################
# Class definitions
#######################################################################

# No classes defined

#######################################################################
# Functions
#######################################################################

def main():
    """Main program
    """
    # Process command line
    p = optparse.OptionParser(usage="%prog OPTIONS gff_file HTSeq_out [ HTSeq_out_2 ... ]",
                              version="%prog "+__version__,
                              description="Annotate HTSeq-count output with data from GFF.")
    p.add_option('-o',action="store",dest="out_file",default=None,
                 help="specify output file name")
    p.add_option('-t','--type',action="store",dest="feature_type",default='exon',
                 help="feature type to process (default 'exon')")
    options,arguments = p.parse_args()
    if len(arguments) < 2:
        p.error("Expected GFF file and at least one HTSeq-count log file")

    # Input files
    gff_file = arguments[0]
    htseq_files = arguments[1:]

    # Feature type being considered
    feature_type = options.feature_type

    # Output file
    if options.out_file:
        annotated_counts_out_file = options.out_file
    else:
        annotated_counts_out_file = os.path.splitext(os.path.basename(gff_file))[0]+"_htseq_counts.txt"
    tables_out_file = os.path.splitext(os.path.basename(annotated_counts_out_file))[0]+"_stats.txt"

    # Process GFF data
    print "Reading data from %s" % gff_file
    gff = GFFcleaner.GFFFile(gff_file)

    # Get a list of unique parent IDs for the feature of interest
    print "Building list of unique parent IDs for '%s' features" % feature_type
    parent_IDs = []
    for data in gff:
        if data['feature'] == feature_type:
            # Get the Parent ID from the attributes
            attributes = GFFcleaner.GFFAttributes(data['attributes'])
            parent = attributes['Parent']
            if parent not in parent_IDs:
                parent_IDs.append(parent)

    # Extract data for features with ID's which match the parents
    print "Creating lookup dictionary for all features with an ID attribute"
    parent_data = {}
    for data in gff:
        # Get the ID from the attributes
        attributes = GFFcleaner.GFFAttributes(data['attributes'])
        try:
            feature_ID = attributes['ID']
        except KeyError:
            # No ID attribute for this feature, skip it
            logging.debug("No ID for line '%s'" % data)
            continue
        if feature_ID in parent_IDs:
            # Extract data for this feature
            try:
                name = attributes['Name']
            except KeyError:
                logging.error("Failed to get name attribute for feature ID %s" % feature_ID)
                sys.exit(1)
            try:
                description = attributes['description']
            except KeyError:
                logging.debug("Failed to get description attribute data for feature ID %s" %
                              feature_ID)
                description = ''
            # Feature type
            parent_feature_type = data['feature']
            # Locus
            locus = "%s:%s-%s" % (data['seqname'],data['start'],data['end'])
            logging.debug("%s\t%s\t%s" % (feature_ID,name,description))
            if feature_ID in parent_data:
                logging.warning("ID '%s' matched multiple times" % parent_ID)
            # Store data
            parent_data[feature_ID] = { 'name': name,
                                        'type': parent_feature_type,
                                        'description': description,
                                        'locus': locus}

    # Process the HTSeq-count files
    print "Processing HTSeq-count files"
    htseq_counts = GFFcleaner.OrderedDictionary()
    htseq_tables = GFFcleaner.OrderedDictionary()
    for htseqfile in htseq_files:
        print "\t%s" % htseqfile
        # Create dictionaries to store data
        htseq_counts[htseqfile] = GFFcleaner.OrderedDictionary()
        htseq_tables[htseqfile] = GFFcleaner.OrderedDictionary()
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

    # Make a list of names found in the count files
    print "Building list of IDs from %s" % htseq_files[0]
    feature_IDs = []
    for feature_ID in htseq_counts[htseq_files[0]]:
        ##print str(feature_ID)
        feature_IDs.append(feature_ID)

    # Create a TabFile for output
    print "Building annotated count file for output"
    header = ['exon_parent',
              'feature_type',
              'locus',
              'name',
              'description']
    for htseqfile in htseq_files:
        header.append(htseqfile)
    annotated_counts = TabFile.TabFile(column_names=header)

    # Combine feature counts and parent feature data
    for feature_ID in feature_IDs:
        # Build the data line
        data = []
        # Add the parent data
        data.append(feature_ID)
        # Look up the name, description etc
        data.append(parent_data[feature_ID]['type'])
        data.append(parent_data[feature_ID]['locus'])
        data.append(parent_data[feature_ID]['name'])
        data.append(parent_data[feature_ID]['description'])
        # Add the counts from each file
        for htseqfile in htseq_files:
            data.append(htseq_counts[htseqfile][feature_ID])
        # Add to the tabfile
        annotated_counts.append(data=data)

    # Write the file
    print "Writing output file %s" % annotated_counts_out_file
    annotated_counts.write(annotated_counts_out_file,include_header=True)

    # Make second file for the trailing table data
    print "Building HTSeq tables file for output"
    header = ['count']
    for htseqfile in htseq_files:
        header.append(htseqfile)
    table_counts = TabFile.TabFile(column_names=header)
    for name in htseq_tables[htseq_files[0]]:
        # Build the data line
        data = [name]
        for htseqfile in htseq_files:
            data.append(htseq_tables[htseqfile][name])
        table_counts.append(data=data)
    print "Writing output file %s" % tables_out_file
    table_counts.write(tables_out_file,include_header=True)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    main()
