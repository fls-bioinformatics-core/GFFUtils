#!/bin/env python
#
#     gffcleaner.py: various functions to clean up GFF files
#     Copyright (C) University of Manchester 2011-2014,2020 Peter Briggs

"""
Utility program which can perform various operations to 'clean' a GFF
file.
"""
#######################################################################
# Module metadata
#######################################################################

from .. import get_version
__version__ = get_version()

#######################################################################
# Import modules
#######################################################################

import os
import sys
import logging
from argparse import ArgumentParser
from ..GFFFile import GFFFile
from ..GFFFile import OrderedDictionary
from ..clean.sgd import GroupByID
from ..clean.sgd import GFFGetDuplicateSGDs
from ..clean.sgd import GFFResolveDuplicateSGDs
from ..clean.sgd import GFFGroupSGDs
from ..clean.sgd import GFFInsertMissingGenes
from ..clean.generic import GFFUpdateAttributes
from ..clean.generic import GFFAddExonIDs
from ..clean.generic import GFFAddIDAttributes
from ..clean.generic import GFFDecodeAttributes
from bcftbx.TabFile import TabFile

# Main program
#
def main():
    """Main program
    """
    # Set up logging format
    logging.basicConfig(format='%(levelname)s: %(message)s')

    p = ArgumentParser(description="Utility to perform various 'cleaning' "
                       "operations on a GFF file")
    p.add_argument('gff_file',metavar="FILE.gff",
                   help="GFF file to operate on")
    p.add_argument('-v','--version',action='version',version=__version__)
    p.add_argument('-o',action='store',dest='output_gff',
                   default=None,
                   help="Name of output GFF file (default is "
                   "'FILE_clean.gff')")
    generic = p.add_argument_group("General cleaning operations")
    generic.add_argument('--prepend',action='store',metavar='STR',
                         dest='prepend_str',default=None,
                         help="String to prepend to seqname in first "
                         "column")
    generic.add_argument('--add-missing-ids',action='store_true',
                         dest='add_missing_ids',default=False,
                         help="For features without an ID attribute, "
                         "construct and insert a generated ID of the "
                         "form '<feature>:<Parent>:<n>' (where n is a "
                         "unique automatically generated number for "
                         "each feature)")
    generic.add_argument('--add-exon-ids',action='store_true',
                         dest='add_exon_ids',default=False,
                         help="For exon features without an ID "
                         "attribute, construct and insert an ID of the "
                         "form 'exon:<Parent>:<n>' (where n is a "
                         "unique automatically generated number for "
                         "each feature)")
    generic.add_argument('--remove-attribute',metavar='ATTR',
                         action='append',dest='rm_attr',
                         help="Remove attribute ATTR from the list of "
                         "attributes for all records in the GFF file. "
                         "Specify this option multiple times to remove "
                         "multiple attributes")
    generic.add_argument('--strict-attributes',action='store_true',
                         dest='strict_attributes',
                         help="Remove attributes that don't conform to "
                         "the KEY=VALUE format specified in the GFF "
                         "standard ")
    generic.add_argument('--no-percent-encoding',action='store_true',
                         dest='no_encoding',default=False,
                         help="Convert encoded attributes to the correct "
                         "characters in the output GFF. !WARNING! this may "
                         "result in a non-cannonical GFF that can be read "
                         "correctly by this or other programs")
    sgdnames = p.add_argument_group("SGD name handling operations")
    sgdnames.add_argument('--report-duplicates',action='store_true',
                          dest='report_duplicates',
                          help="Report features with duplicated SGD "
                          "names (i.e. identical values for the 'SGD' "
                          "attribute) and write list to "
                          "'FILE_duplicates.gff' with line numbers, "
                          "chromosome, start coordinate and strand")
    sgdnames.add_argument('--resolve-duplicates',action='store',
                          dest='mapping_file',default=None,
                          help="Resolve presence of features with "
                          "duplicate SGD names, by selecting a single "
                          "feature determined by matching against "
                          "'best' genes in the supplied MAPPING_FILE "
                          "(tab-delimited file with columns "
                          "'gene_name,chromosome,start,end,strand'). "
                          "If a best match can be found then it is kept "
                          "while remaining non-matching duplicate "
                          "features are removed and written to "
                          "'FILE_discarded.gff'; otherwise the "
                          "duplication is unresolved and all features "
                          "are kept (use --discard-unresolved to remove "
                          "these)")
    sgdnames.add_argument('--discard-unresolved',action='store_true',
                          dest='discard_unresolved',
                          help="Remove all features with duplicated "
                          "SGD names which cannot be resolved by "
                          "--resolve-duplicates; the discarded features "
                          "will be written to 'FILE_unresolved.gff'")
    sgdnames.add_argument('--insert-missing',action='store',dest='gene_file',
                          default=None,
                          help="Insert genes from GENE_FILE with names "
                          "that don't appear in any SGD attributes of "
                          "features in the input GFF. GENE_FILE should "
                          "be a tab-delimited file with columns "
                          "'gene_name,chromosome,start,end,strand'. (If "
                          "no GENE_FILE is supplied then the MAPPING_FILE "
                          "supplied via --resolve-duplicates will be "
                          "used)")
    sgdops = p.add_argument_group("SGD-specific cleaning operations")
    sgdops.add_argument('--clean-exclude-attributes',action='store_true',
                        dest='do_clean_exclude_attributes',
                        help="Remove 'kaks', 'kaks2' and 'ncbi' attributes, "
                        "if present. (To remove arbitrary attributes, use "
                        "--remove-attribute)")
    sgdops.add_argument('--clean-group-sgds',action='store_true',
                        dest='do_clean_group_sgds',
                        help="Group features with the same SGD by adding "
                        "unique numbers to the 'ID' attributes; IDs will "
                        "have the form 'CDS:<SGD>:<n>' (where n is a unique "
                        "number for a given SGD)")
    sgdops.add_argument('--clean-replace-attributes',action='store_true',
                        dest='do_clean_replace_attributes',
                        help="Replace 'ID', 'Gene', 'Parent' and 'Name' "
                        "attributes with the value of the SGD attribute, "
                        "if present")
    sgdops.add_argument('--clean-score',action='store_true',
                        dest='do_clean_score',
                        help="Replace 'Anc_*' and blanks in 'score' field "
                        "with zeroes")
    sgdops.add_argument('--clean',action='store_true',dest='do_clean',
                        help="Perform all the 'cleaning' manipulations on "
                        "the input data (equivalent to specifying all of "
                        "--clean-score, --clean-replace-attributes, "
                        "--clean-exclude-attributes and --clean-group-sgds)")
    advanced = p.add_argument_group("Advanced options")
    advanced.add_argument('--debug',action='store_true',dest='debug',
                          help="Print debugging information")

    # Process the command line
    args = p.parse_args()

    # Check for debugging
    if args.debug:
        # Turn on debugging output
        logging.getLogger().setLevel(logging.DEBUG)

    # Input file
    infile = args.gff_file
    if not os.path.exists(infile):
        p.error("Input file '%s' not found" % infile)

    # Report version
    p.print_version()

    # Set flags based on command line

    # String to prepend to first column
    prepend_str = args.prepend_str
    # Cleaning options
    if args.do_clean:
        # Update values in the "score" column
        clean_score = True
        # Clean up the "attributes" column
        clean_replace_attributes = True
        clean_exclude_attributes = True
        # Set ID field in "attributes" to group lines with matching SGDs
        group_SGDs = True
    else:
        # Set options based on user input
        clean_score = args.do_clean_score
        clean_replace_attributes = args.do_clean_replace_attributes
        clean_exclude_attributes = args.do_clean_exclude_attributes
        group_SGDs = args.do_clean_group_sgds
    # Report duplicate names
    report_duplicates = args.report_duplicates
    # Resolve duplicated genes using CDS file
    if args.mapping_file is not None:
        resolve_duplicates = True
        cdsfile = args.mapping_file
    else:
        resolve_duplicates = False
        cdsfile = None
    # Discard unresolved duplicates
    discard_unresolved = args.discard_unresolved
    # Insert missing genes
    if args.gene_file is not None:
        insert_missing = True
        if args.gene_file:
            genefile = args.gene_file
        else:
            genefile = cdsfile
    else:
        insert_missing = False
        genefile = None
    # Add an artificial exon ID attribute
    add_exon_ids = args.add_exon_ids
    # Add generated ID attributes
    add_missing_ids = args.add_missing_ids
    # Suppress encoding of attributes on output
    no_attribute_encoding = args.no_encoding
    # Remove attributes that don't conform to KEY=VALUE format
    strict_attributes = args.strict_attributes

    # Name for output files
    if not args.output_gff:
        outbase = os.path.splitext(os.path.basename(infile))[0]
        outfile = outbase+'_clean.gff'
    else:
        outbase = os.path.splitext(os.path.basename(args.output_gff))[0]
        outfile = args.output_gff
    print("Input : %s" % infile)
    print("Output: %s" % outfile)
    dupfile = outbase+'_duplicates.txt'
    delfile = outbase+'_discarded.gff'
    unresfile = outbase+'_unresolved.gff'

    # Read in data from file
    gff_data = GFFFile(infile)

    # Prepend string to seqname column
    if prepend_str is not None:
        print("Prepending '%s' to values in 'seqname' column" % prepend_str)
        for data in gff_data:
            data['seqname'] = prepend_str+str(data['seqname'])

    # Check/clean score column values
    if clean_score:
        print("Replacing 'Anc_*' and blanks with '0's in 'score' column")
        score_unexpected_values = set()
        for data in gff_data:
            try:
                # Numerical value
                score = float(data['score'])
                if score != 0:
                    score_unexpected_values.add(data['score'])
            except ValueError:
                # String value
                if data['score'].startswith('Anc_') or \
                   data['score'].strip() == '':
                    # Replace "Anc_*" or blank values in "score"
                    # column with zero
                    data['score'] = '0'
                else:
                    score_unexpected_values.add(data['score'])
        # Report unexpected values
        score_unexpected_values = sorted(list(score_unexpected_values))
        n = len(score_unexpected_values)
        if n > 0:
            logging.warning("%d 'score' values that are not '', 0 or 'Anc_*'" % n)
            logging.warning("Other values: %s" %
                            ', '.join([str(x)
                                       for x in score_unexpected_values]))

    # Clean up the data in "attributes" column: replace keys
    if clean_replace_attributes:
        # Initialise mapping of keys from input to output in "attributes" column
        # where new values are required etc
        attributes_key_map = OrderedDictionary()
        attributes_key_map['ID'] = 'SGD'
        attributes_key_map['Gene'] = 'SGD'
        attributes_key_map['Parent'] = 'SGD'
        attributes_key_map['Name'] = 'SGD'
        attributes_dont_replace_with_empty_data = True
        print("Cleaning up attributes: replacing keys:")
        for key in attributes_key_map.keys():
            print("\t%s -> %s" % (key,attributes_key_map[key]))
        if attributes_dont_replace_with_empty_data:
            print("(Replacement will be skipped if new data is missing/blank)")
        GFFUpdateAttributes(gff_data,attributes_key_map,[],
                            attributes_dont_replace_with_empty_data)

    # Clean up the data in "attributes" column: exclude keys
    if clean_exclude_attributes:
        # List of keys to exclude
        attributes_exclude_keys = ['kaks','kaks2','ncbi']
        print("Excluding keys:")
        for key in attributes_exclude_keys:
            print("\t%s" % key)
        GFFUpdateAttributes(gff_data,{},attributes_exclude_keys,True)

    # Set the IDs for consecutive lines with matching SGD names, to
    # indicate that they're in the same gene
    if group_SGDs:
        print("Grouping SGDs by setting ID's for consecutive lines "
              "with the same SGD values")
        GFFGroupSGDs(gff_data)

    # Find duplicates in input file
    if report_duplicates or resolve_duplicates:
        duplicate_sgds = GFFGetDuplicateSGDs(gff_data)
                
    if report_duplicates:
        # Write to duplicates file
        print("Writing duplicate SGD names to %s" % dupfile)
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
            for group in GroupByID(duplicate_sgds[sgd]):
                if len(group) > 1: ngroups += 1
        if ndup == 0:
            fd.write("No duplicate SGDs\n")
        fd.close()
        print("%d duplicates found (of which %d are trivial)" %
              (ndup,ngroups))

    if resolve_duplicates:
        print("Resolving duplicate SGDs using data from %s" % cdsfile)
        print("Discarded genes will be written to %s" % delfile)
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
            print("No mapping genes with same SGDs found in %s:" % cdsfile)
            for sgd in unresolved_sgds_no_mapping_genes:
                print("\t%s" % sgd)
            print("")
        if len(unresolved_sgds_no_mapping_genes_after_filter) > 0:
            print("No mapping genes with same chromosome and/or strand:")
            for sgd in unresolved_sgds_no_mapping_genes_after_filter:
                print("\t%s" % sgd)
            print("")
        if len(unresolved_sgds_no_overlaps) > 0:
            print("No mapping genes with overlaps:")
            for sgd in unresolved_sgds_no_overlaps:
                print("\t%s" % sgd)
            print("")
        if len(unresolved_sgds_multiple_matches) > 0:
            print("Multiple matching mapping genes:")
            for sgd in unresolved_sgds_multiple_matches:
                print("\t%s" % sgd)
            print("")
        # Summary counts for each case
        print("Total number of duplicated indexes   : %d" %
              len(duplicate_sgds.keys()))
        print("Number of resolved duplicate SGDs    : %d" %
              len(resolved_sgds))
        print("Unresolved duplicates:")
        print("* No mapping genes with same SGD     : %d" %
              len(unresolved_sgds_no_mapping_genes))
        print("* No mapping genes with same chr/str : %d" %
              len(unresolved_sgds_no_mapping_genes_after_filter))
        print("* No mapping genes with overlap      : %d" %
              len(unresolved_sgds_no_overlaps))
        print("* Multiple mapping genes match       : %d" %
              len(unresolved_sgds_multiple_matches))

        # Remove discarded duplicates from the data
        print("Removing discarded duplicates and writing to %s" % delfile)
        fd = open(delfile,'w')
        for discard_data in discard:
            try:
                ip = gff_data.indexByLineNumber(discard_data.lineno())
                del(gff_data[ip])
                fd.write("%s\n" % discard_data)
            except IndexError:
                logging.warning("Failed to delete line %d: not found" %
                                discard_data.lineno())
        fd.close()

        # Remove unresolved duplicates if requested
        if discard_unresolved:
            print("Removing unresolved duplicates and writing to %s" %
                  unresfile)
            # Get list of unresolved SGDs
            all_unresolved = result['unresolved_sgds']
            # Get list of unresolved duplicates
            unresolved = []
            for data in gff_data:
                attributes = data['attributes']
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
        print("Inserting unmatched genes from %s" % genefile)
        # Get gene data from CDS file
        # Format is tab-delimited, each line has:
        # orf      chr      start     end      strand
        mapping = TabFile(genefile,column_names=('name','chr','start','end','strand'))
        n_genes_before_insert = len(gff_data)
        gff_data = GFFInsertMissingGenes(gff_data,mapping)
        print("Inserted %d missing genes" %
              (len(gff_data) - n_genes_before_insert))

    # Construct and insert ID for exons
    if add_exon_ids:
        print("Inserting artificial IDs for exon records")
        gff_data = GFFAddExonIDs(gff_data)

    # Construct and insert missing ID attributes
    if add_missing_ids:
        print("Inserting generated IDs for records where IDs are missing")
        gff_data = GFFAddIDAttributes(gff_data)

    # Strip attributes requested for removal
    if args.rm_attr:
        print("Removing the following attributes from all records:")
        for attr in args.rm_attr:
            print("\t* %s" % attr)
        GFFUpdateAttributes(gff_data,exclude_keys=args.rm_attr)

    # Remove attributes that don't conform to KEY=VALUE format
    if strict_attributes:
        print("Removing attributes that don't conform to KEY=VALUE format")
        GFFUpdateAttributes(gff_data,exclude_nokeys=True)

    # Suppress percent encoding of attributes
    if no_attribute_encoding:
        print("Converting encoded special characters in attribute data to "
              "non-encoded form")
        logging.warning("!!! Special characters will not be correctly "
                        "encoded in the output  !!!")
        logging.warning("!!! The resulting GFF may not be readable by this "
                        "or other programs !!!")
        gff_data = GFFDecodeAttributes(gff_data)

    # Write to output file
    print("Writing output file %s" % outfile)
    gff_data.write(outfile)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s: %(message)s")
    main()


