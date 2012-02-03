#!/usr/bin/env python
#
#     gff_extract_mRNA_to_gene.py: match gene data to mRNA from GFF file
#     Copyright (C) University of Manchester 2012 Peter Briggs, Leo Zeef
#
#######################################################################
#
# gff_extract_mRNA_to_gene.py
#
#######################################################################

"""gff_extract_mRNA_to_gene.py

Utility to match 'mRNA' features with 'gene' features in a GFF, by
finding genes which have the same 'ID' attribute as the mRNA's 'Parent'
attribute. Outputs a list of mRNA IDs with the matching gene ID, name
and description.

Uses the GFFFile and supporting classes to read in the GFF file.
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
    p = optparse.OptionParser(usage="%prog OPTIONS gff_file",
                              version="%prog "+__version__,
                              description="Match mRNA IDs in a GFF file with their 'parent' "
                              "genes: reads in GFF file and for each 'mRNA' feature, identifies "
                              "matching genes where the gene 'ID' attribute is the same as the "
                              "mRNA 'Parent' attribute. Outputs a tab-delimited file with each "
                              "line consisting of mRNA ID, gene ID, gene name and gene "
                              "description. "
                              "Output file name can be specified with the -o option, otherwise "
                              "it will be the input file name with '_annotated' appended.")
    p.add_option('-o',action="store",dest="out_file",default=None,
                 help="specify output file name")
    options,arguments = p.parse_args()
    if len(arguments) != 1:
        p.error("Input GFF file expected")

    # Input and output file
    gff_file = arguments[0]
    if options.out_file:
        out_file = options.out_file
    else:
        out_file = os.path.splitext(os.path.basename(gff_file))[0]+"_mRNA_to_gene.txt"
        print "Output file: %s" % out_file

    # Read in the GFF
    gff = GFFcleaner.GFFFile(gff_file)

    # Assemble list of gene IDs from mRNA lines, and store data for
    # genes (ID, Name and Description)
    mRNA_to_genes = {}
    gene_data = {}
    for line in gff:
        if line['feature'] == "mRNA":
            # Get the mRNA ID and the Parent (= gene ID) from the attributes
            attributes = GFFcleaner.GFFAttributes(line['attributes'])
            mRNA_ID = attributes['ID']
            gene_ID = attributes['Parent']
            logging.debug("mRNA_ID = %s\tgene_ID = %s" % (mRNA_ID,gene_ID))
            mRNA_to_genes[mRNA_ID] = gene_ID
        elif line['feature'] == "gene":
            # Get the gene ID, Name and Description
            attributes = GFFcleaner.GFFAttributes(line['attributes'])
            gene_ID = attributes['ID']
            try:
                name = attributes['Name']
            except KeyError:
                logging.error("Failed to get name attribute for gene ID %s" % gene_ID)
                sys.exit(1)
            try:
                description = attributes['description']
            except KeyError:
                logging.debug("Failed to get description attribute data for gene ID %s" % gene_ID)
                description = ''
            logging.debug("%s\t%s\t%s" % (gene_ID,name,description))
            if gene_ID in gene_data:
                logging.warning("gene ID '%s' matched multiple times" % gene_ID)
            gene_data[gene_ID] = { 'name': name, 'description': description }

    # Combine mRNA and gene data for output
    fo = open(out_file,'w')
    for mRNA_ID in mRNA_to_genes:
        gene_ID = mRNA_to_genes[mRNA_ID]
        # mRNA ID and gene ID (i.e. parent)
        data = [mRNA_ID,gene_ID]
        # Append gene data
        data.append(gene_data[gene_ID]['name'])
        data.append(gene_data[gene_ID]['description'])
        fo.write('\t'.join(data)+'\n')
    fo.close()

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    main()
