#!/usr/bin/env python

import optparse
import GFFcleaner
import sys

__version__ = "0.0.1"

if __name__ == "__main__":
    # Process command line
    p = optparse.OptionParser(usage="%prog OPTIONS gff_file",
                              version="%prog "+__version__,
                              description="")
    options,arguments = p.parse_args()
    if len(arguments) != 1:
        p.error("Input GFF file expected")

    # Input file
    gff_file = arguments[0]

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
            ##print "mRNA_ID = %s\tgene_ID = %s" % (mRNA_ID,gene_ID)
            mRNA_to_genes[mRNA_ID] = gene_ID
        elif line['feature'] == "gene":
            # Get the gene ID, Name and Description
            attributes = GFFcleaner.GFFAttributes(line['attributes'])
            gene_ID = attributes['ID']
            try:
                name = attributes['Name']
            except KeyError:
                print "Failed to get name attribute for gene ID %s" % gene_ID
                sys.exit(1)
            try:
                description = attributes['description']
            except KeyError:
                ##print "Failed to get description attribute data for gene ID %s" % gene_ID
                description = ''
            ##print "%s\t%s\t%s" % (gene_ID,name,description)
            gene_data[gene_ID] = { 'name': name, 'description': description }

    # Combine mRNA and gene data for output
    for mRNA_ID in mRNA_to_genes:
        gene_ID = mRNA_to_genes[mRNA_ID]
        # mRNA ID and gene ID (i.e. parent)
        data = [mRNA_ID,gene_ID]
        # Append gene data
        data.append(gene_data[gene_ID]['name'])
        data.append(gene_data[gene_ID]['description'])
        print '\t'.join(data)
        
        
            

            

