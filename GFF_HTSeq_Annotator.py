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

Usage
-----

First Run the htseq-count program to generate counts using e.g.

   htseq-count -q -i ID -t exon <sam> <gff> > htseq-counts

Then run the annotator:

   GFF_HTSeq_Annotator.py -t exon <gff> htseq-counts

Multiple parents
----------------

It's possible for features in a GFF file to have multiple parents.

In this case the output from htseq-count will reproduce the 'Parent'
attribute verbatim, e.g. AF2312,AB2812,abc-3.

GFF_HTSeq_Annotator will 

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
                              description="Annotate htseq-count output with data from GFF. "
                              "Generate input HTSeq_out files using the htseq-count program "
                              "e.g. 'htseq-count -q -t exon -i Parent gff_file sam_file'. "
                              "The annotator looks up the parent genes of each feature and "
                              "outputs this information against the htseq-counts (in"
                              "<gff_file>_htseq_counts.txt) plus the totals assigned, not "
                              "counted etc (in <gff_file>_htseq_counts_stats.txt).")
    p.add_option('-o',action="store",dest="out_file",default=None,
                 help="specify output file name")
    p.add_option('-t','--type',action="store",dest="feature_type",default='exon',
                 help="feature type to process (default 'exon'; should be the same as the "
                 "feature type used in htseq-count runs)")
    options,arguments = p.parse_args()
    if len(arguments) < 2:
        p.error("Expected GFF file and at least one HTSeq-count log file")

    # Input GFF file
    gff_file = arguments[0]
    if not os.path.exists(gff_file):
        p.error("Input GFF file %s not found" % gff_file)

    # Check for wildcards in HTSeq file names, to emulate linux shell globbing
    # on platforms such as Windows which don't have this built in
    htseq_files = []
    for arg in arguments[1:]:
        for filen in glob.iglob(arg):
            if not os.path.exists(filen):
                p.error("File '%s' not found" % filen)
            htseq_files.append(filen)
    if not htseq_files:
        p.error("No input HTSeq-count files found")

    # Feature type being considered
    feature_type = options.feature_type

    # Output files
    if options.out_file:
        annotated_counts_out_file = options.out_file
    else:
        annotated_counts_out_file = os.path.splitext(os.path.basename(gff_file))[0]+\
            "_htseq_counts.txt"
    tables_out_file = \
        os.path.splitext(os.path.basename(annotated_counts_out_file))[0]+\
        "_stats"+os.path.splitext(annotated_counts_out_file)[1]

    # Process GFF data
    print "Reading data from %s" % gff_file
    gff = GFFFile.GFFFile(gff_file)

    # Overview of the method:
    #
    # HTSeq-count output has lists of exon parent IDs plus counts
    #
    # GFF hierarchy has exon -> [intermediate] -> gene
    # where intermediate could be mRNA, pseudogene etc
    #
    # For all the exon parent IDs from HTSeq-count, locate the intermediates and then
    # for each intermediate locate the genes.

    # Get a list of unique parent IDs for the feature of interest (=exon by default)
    print "Building list of unique parent IDs for '%s' features" % feature_type
    # List of exon parent IDs
    exon_parent_IDs = []
    for data in gff:
        if data['feature'] == feature_type:
            # Get the Parent ID from the attributes
            attributes = GFFFile.GFFAttributes(data['attributes'])
            parent = attributes['Parent']
            # Check for multiple parents - split on comma
            if len(parent.split(',')) > 1:
                # Issue a warning but continue for now
                logging.warning("Multiple parents found on line %d: %s" % (data.lineno(),
                                                                           parent))
            # Store the ID
            if parent not in exon_parent_IDs:
                exon_parent_IDs.append(parent)

    # Get the actual data for each exon parent
    # Also build a list of *their* parents (i.e. genes)
    print "Collecting data for %s parents & building lookup list of their parent gene IDs" % \
        feature_type
    # Data about exon parents (e.g. type, the ID of its parent)
    exon_parent_data = {}
    # Link IDs of genes to the exon parents that they (the genes) are parents of
    parent_gene_IDs = {}
    for data in gff:
        attributes = GFFFile.GFFAttributes(data['attributes'])
        try:
            feature_ID = attributes['ID']
            if feature_ID in exon_parent_IDs:
                # This is a parent of an exon - store its data
                feature_parent_ID = attributes['Parent']
                # Check for multiple parents
                if len(feature_parent_ID.split(',')) > 1:
                    # Issue a warning
                    logging.warning("Multiple parents detected for %s ID %s (line %d,parents %s)" %
                                    (data['feature'],feature_ID,data.lineno(),feature_parent_ID))
                exon_parent_data[feature_ID] = { 'Parent': feature_parent_ID,
                                                 'type' : data['feature'] }
                # Also store lookup from gene IDs to exon parents
                parent_gene_IDs[feature_parent_ID] = feature_ID
        except KeyError:
            # No ID attribute for this feature, skip it
            logging.debug("No ID for line '%s'" % data)

    # Extract data for features with ID's which match the exon parent's parents (i.e. genes)
    print "Collecting data for parent genes"
    # Data about the genes (e.g. name, locus, description)
    parent_genes = {}
    for data in gff:
        # Get the ID from the attributes
        attributes = GFFFile.GFFAttributes(data['attributes'])
        try:
            feature_ID = attributes['ID']
        except KeyError:
            # No ID attribute for this feature, skip it
            logging.debug("No ID for line '%s'" % data)
            continue
        if feature_ID in parent_gene_IDs:
            # This is one of the exon parent's parent genes: extract and store the data
            try:
                name = attributes['Name']
            except KeyError:
                logging.error("Failed to get name attribute for feature ID %s" % feature_ID)
                sys.exit(1)
            # Check it's a gene
            if data['feature'] != 'gene':
                logging.error("Non-gene parent! %s" % feature_ID)
                sys.exit(1)
            # Build description text
            # This will be all the attribute data from the 'description' attribute
            # onwards
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
            # Locus: chromosome plus start and end data
            locus = "%s:%s-%s" % (data['seqname'],data['start'],data['end'])
            logging.debug("%s\t%s\t%s" % (feature_ID,name,description))
            # Gene length
            gene_length = data['end'] - data['start']
            # Check if it's unique
            if feature_ID in parent_genes:
                logging.warning("ID '%s' matched multiple times" % feature_ID)
            # Store data
            parent_genes[feature_ID] = { 'gene_name': name,
                                         'gene_locus': locus,
                                         'description': description,
                                         'chr': data['seqname'],
                                         'start': data['start'],
                                         'end': data['end'],
                                         'gene_length': gene_length }

    # Process the HTSeq-count files
    print "Processing HTSeq-count files"
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
        ##print str(feature_ID)
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
                                                     'gene_length',
                                                     'locus',
                                                     'description'])
    for htseqfile in htseq_files:
        annotated_counts.appendColumn(htseqfile)

    # Combine feature counts and parent feature data
    for exon_parent_ID in feature_IDs:
        # Check that we have data for the parent referenced in the HTSeq_count files
        if exon_parent_ID in exon_parent_data:
            # Add the exon parent data
            data.append(exon_parent_ID)
            exon_parent_gene = exon_parent_data[exon_parent_ID]['Parent']
            # Look up the name, description etc
            exon_parent_type        = exon_parent_data[exon_parent_ID]['type']
            exon_parent_gene_name   = parent_genes[exon_parent_gene]['gene_name']
            exon_parent_chr         = parent_genes[exon_parent_gene]['chr']
            exon_parent_start       = parent_genes[exon_parent_gene]['start']
            exon_parent_end         = parent_genes[exon_parent_gene]['end']
            exon_parent_gene_length = parent_genes[exon_parent_gene]['gene_length']
            exon_parent_gene_locus  = parent_genes[exon_parent_gene]['gene_locus']
            exon_parent_description = parent_genes[exon_parent_gene]['description']
        else:
            # No data found for the exon parent 
            logging.error("No data found in GFF for %s parent '%s' referenced in "
                          "HTSeq_count files" % (feature_type,exon_parent_ID))
            exon_parent_type        = ''
            exon_parent_gene_name   = ''
            exon_parent_chr         = ''
            exon_parent_start       = ''
            exon_parent_end         = ''
            exon_parent_gene_length = ''
            exon_parent_gene_locus  = ''
            exon_parent_description = ''
        # Build the data line
        data = [exon_parent_ID,
                exon_parent_type,
                exon_parent_gene_name,
                exon_parent_chr,
                exon_parent_start,
                exon_parent_end,
                exon_parent_gene_length,
                exon_parent_gene_locus,
                exon_parent_description]
        # Add the counts from each file
        for htseqfile in htseq_files:
            data.append(htseq_counts[htseqfile][exon_parent_ID])
        # Add to the tabfile
        annotated_counts.append(data=data)

    # Write the file
    print "Writing output file %s" % annotated_counts_out_file
    annotated_counts.write(annotated_counts_out_file,include_header=True,no_hash=True)

    # Make second file for the trailing table data
    print "Building HTSeq tables file for output"
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
