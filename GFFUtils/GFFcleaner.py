#!/bin/env python
#
#     GFFcleaner.py: various functions to clean up GFF files
#     Copyright (C) University of Manchester 2011-2014 Peter Briggs
#
######################################################################
#
# GFFcleaner.py
#
#######################################################################

"""GFFcleaner

Utility program which can perform various operations to 'clean' a GFF
file.
"""
#######################################################################
# Module metadata
#######################################################################

from . import get_version
__version__ = get_version()

#######################################################################
# Import modules
#######################################################################

import os,sys
import logging
import optparse

try:
    from TabFile import TabFile
except ImportError:
    from bcftbx.TabFile import TabFile
from GFFFile import GFFFile,GFFAttributes,GFFID,OrderedDictionary

#######################################################################
# Classes
#######################################################################

# No classes defined

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
        this_index = GFFID(data['attributes']['ID']).index
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
            logging.debug("\t%s" % GFFID(gene['attributes']['ID']))
    return subsets

def GFFUpdateAttributes(gff_data,update_keys={},exclude_keys=[],no_empty_values=True,
                        exclude_nokeys=False):
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
      exclude_nokeys: if True then any 'nokeys' attributes will be removed
    """
    for data in gff_data:
        # Process the attributes data
        attributes = data['attributes']
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
                    logging.warning("Cannot update value of attribute '%s': "
                                    "replacement attribute '%s' not found" %
                                    (key,lookup_key))
            except KeyError:
                # No mapping found for key, ignore
                pass
        # Remove 'nokeys' data
        if exclude_nokeys:
            del(attributes.nokeys()[:])
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
        attributes = data['attributes']
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
                attr = duplicate['attributes']
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
        attributes = data['attributes']
        # Get the SGD value
        try:
            sgd = attributes['SGD']
        except KeyError:
            # SGD not in the attributes, treat as blank
            sgd = ''
        if sgd != '':
            # Check the ID
            idx = GFFID(attributes['ID'])
            if idx.code != 'CDS':
                # Set the CDS prefix and index and update ID attribute
                idx.code = 'CDS'
                idx.index = 1
                attributes['ID'] = str(idx)
            ln = data.lineno()
            # Loop over next 5 data lines after this looking for matching SGD
            for data0 in gff_data[next_ln:next_ln+5]:
                attr0 = data0['attributes']
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
        attributes = data['attributes']
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
            attributes = missing['attributes']
            attributes['ID'] = 'CDS:%s:1' % sgd
            attributes['SGD'] = sgd
            attributes['Gene'] = sgd
            attributes['Parent'] = sgd
    # Finished inserting missing genes
    return gff_data

def GFFAddExonIDs(gff_data):
    """Construct and insert a ID attribute for exons

    For each exon in the input GFF data, insert an ID attribute of
    the form:

    ID=exon:<Parent>:<n>

    where <Parent> is the name specified in the Parent attribute and
    <n> is a numerical string that is unique across all exons in the
    GFF data.

    Note that if an exon already has an ID attribute then it will be
    overwritten with the constructed ID string.

    Arguments:
      gff_data: a GFFFile object containing the GFF file data (which
        is modified in place)

    Returns:
      The modified GFFFile object.
    """
    count = 0
    for record in gff_data:
        if record['feature'] == 'exon':
            attributes = record['attributes']
            if 'Parent' not in attributes:
                logging.warning("No 'Parent' attribute")
            else:
                count += 1
                exon_ID = "exon:%s:%08d" % (attributes['Parent'],count)
                if 'ID' not in attributes:
                    attributes.insert(0,'ID',exon_ID)
                else:
                    attributes['ID'] = exon_ID
    return gff_data

def GFFAddIDAttributes(gff_data):
    """Construct and insert a ID attribute for all features without one

    For each feature in the input GFF data that doesn't already have an
    an ID attribute, insert one of the form:

    ID=<feature>:<Parent>:<n>

    where <feature> is the feature type (e.g. 'exon', 'CDS' etc),
    <Parent> is the name specified in the Parent attribute and
    <n> is a numerical string that is unique across all exons in the
    GFF data.

    Arguments:
      gff_data: a GFFFile object containing the GFF file data (which
        is modified in place)

    Returns:
      The modified GFFFile object.
    """
    count = 0
    for record in gff_data:
        attributes = record['attributes']
        if 'ID' not in attributes:
            # Add an ID
            count += 1
            if 'Parent' not in attributes:
                logging.warning("No 'Parent' attribute")
                feature_ID = "%s:NOPARENT:%08d" % (record['feature'],
                                                   count)
            else:
                feature_ID = "%s:%s:%08d" % (record['feature'],
                                             attributes['Parent'],
                                             count)
                attributes.insert(0,'ID',feature_ID)
    return gff_data

def GFFDecodeAttributes(gff_data):
    """Update GFF records removing percent encoding of special characters

    The GFF format specifies that certain special characters in the
    "attribute" field must be escaped using "percent encoding" (i.e.
    HTML-style encoding).

    This function converts special characters back to their original
    form, essentially 'decoding' the attribute values so that they are
    more "human-readable". However the resulting attributes may be non-
    cannonical and thus cause problems when fed into other GFF reading
    programs that parse the attribute field. Caveat emptor.

    Arguments:
      gff_data: a GFFFile object containing the GFF file data (which
        is modified in place)

    Returns:
      The modified GFFFile object.
    """
    for record in gff_data:
        attributes = record['attributes']
        attributes.encode(False)
    return gff_data

# Main program
#
def main():
    """Main program
    """
    # Set up logging format
    logging.basicConfig(format='%(levelname)s: %(message)s')

    p = optparse.OptionParser(usage="%prog [options] <file>.gff",
                              version="%prog "+__version__,
                              description=
                              "Utility to perform various 'cleaning' operations on a GFF file.")
    p.add_option('-o',action='store',dest='output_gff',default=None,
                 help="Name of output GFF file (default is '<file>_clean.gff')")
    p.add_option('--prepend',action='store',dest='prepend_str',default=None,
                 help="String to prepend to seqname in first column")
    p.add_option('--clean',action='store_true',dest='do_clean',
                 help="Perform all the 'cleaning' manipulations on the input data (equivalent "
                 "to specifying all of --clean-score, --clean-replace-attributes, "
                 "--clean-exclude-attributes and --clean-group-sgds)")
    p.add_option('--clean-score',action='store_true',dest='do_clean_score',
                 help="Replace 'Anc_*' and blanks in 'score' field with zeroes")
    p.add_option('--clean-replace-attributes',action='store_true',
                 dest='do_clean_replace_attributes',
                 help="Replace 'ID', 'Gene', 'Parent' and 'Name' attributes with the value "
                 "of the SGD attribute, if present")
    p.add_option('--clean-exclude-attributes',action='store_true',
                 dest='do_clean_exclude_attributes',
                 help="Remove the 'kaks', 'kaks2' and 'ncbi' attributes (to remove "
                 "arbitrary attributes, see the --remove-attribute=... option)")
    p.add_option('--clean-group-sgds',action='store_true',dest='do_clean_group_sgds',
                 help="Group features with the same SGD by adding unique numbers to the 'ID' "
                 "attributes; IDs will have the form 'CDS:<SGD>:<n>' (where n is a unique "
                 "number for a given SGD)")
    p.add_option('--report-duplicates',action='store_true',dest='report_duplicates',
                 help="Report duplicate SGD names and write list to <file>_duplicates.gff "
                 "with line numbers, chromosome, start coordinate and strand.")
    p.add_option('--resolve-duplicates',action='store',dest='mapping_file',default=None,
                 help="Resolve duplicate SGDs by matching against 'best' genes in the supplied "
                 "mapping file; other non-matching genes are discarded and written to "
                 "<file>_discarded.gff.")
    p.add_option('--discard-unresolved',action='store_true',dest='discard_unresolved',
                 help="Discard any unresolved duplicates, which are written to "
                 "<file>_unresolved.gff.")
    p.add_option('--insert-missing',action='store',dest='gene_file',default=None,
                 help="Insert genes from gene file with SGD names that don't appear in the "
                 "input GFF. If GENE_FILE is blank ('='s must still be present) then the mapping "
                 "file supplied with the --resolve-duplicates option will be used instead.")
    p.add_option('--add-exon-ids',action='store_true',dest='add_exon_ids',default=False,
                 help="For exon features without an ID attribute, construct and insert an "
                 "ID of the form 'exon:<Parent>:<n>' (where n is a unique number).")
    p.add_option('--add-missing-ids',action='store_true',dest='add_missing_ids',default=False,
                 help="For features without an ID attribute, construct and insert a "
                 "generated ID of the form '<feature>:<Parent>:<n>' (where n is a unique "
                 "number).")
    p.add_option('--no-percent-encoding',action='store_true',dest='no_encoding',default=False,
                 help="Convert encoded attributes to the correct characters in "
                 "the output GFF. WARNING this may result in a non-cannonical GFF that can't "
                 "be read correctly by this or other programs.")
    p.add_option('--remove-attribute',action='append',dest='rm_attr',
                 help="Remove attribute RM_ATTR from the list of attributes for all records "
                 "in the GFF file (can be specified multiple times)")
    p.add_option('--strict-attributes',action='store_true',dest='strict_attributes',
                 help="Remove attributes that don't conform to the KEY=VALUE format")
    p.add_option('--debug',action='store_true',dest='debug',
                 help="Print debugging information")

    # Process the command line
    options,arguments = p.parse_args()

    # Check for debugging
    if options.debug:
        # Turn on debugging output
        logging.getLogger().setLevel(logging.DEBUG)

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
        clean_replace_attributes = True
        clean_exclude_attributes = True
        # Set ID field in "attributes" to group lines with matching SGDs
        group_SGDs = True
    else:
        # Set options based on user input
        clean_score = options.do_clean_score
        clean_replace_attributes = options.do_clean_replace_attributes
        clean_exclude_attributes = options.do_clean_exclude_attributes
        group_SGDs = options.do_clean_group_sgds
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
    # Add an artificial exon ID attribute
    add_exon_ids = options.add_exon_ids
    # Add generated ID attributes
    add_missing_ids = options.add_missing_ids
    # Suppress encoding of attributes on output
    no_attribute_encoding = options.no_encoding
    # Remove attributes that don't conform to KEY=VALUE format
    strict_attributes = options.strict_attributes

    # Name for output files
    ##outext = os.path.splitext(os.path.basename(infile))[1]
    if not options.output_gff:
        outbase = os.path.splitext(os.path.basename(infile))[0]
        outfile = outbase+'_clean.gff'
    else:
        outbase = os.path.splitext(os.path.basename(options.output_gff))[0]
        outfile = options.output_gff
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
        print "Cleaning up attributes: replacing keys:"
        for key in attributes_key_map.keys():
            print "\t%s -> %s" % (key,attributes_key_map[key])
        if attributes_dont_replace_with_empty_data:
            print "(Replacement will be skipped if new data is missing/blank)"
        GFFUpdateAttributes(gff_data,attributes_key_map,[],
                            attributes_dont_replace_with_empty_data)

    # Clean up the data in "attributes" column: exclude keys
    if clean_exclude_attributes:
        # List of keys to exclude
        attributes_exclude_keys = ['kaks','kaks2','ncbi']
        print "Excluding keys:"
        for key in attributes_exclude_keys:
            print "\t%s" % key
        GFFUpdateAttributes(gff_data,{},attributes_exclude_keys,True)

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
        print "Inserting unmatched genes from %s" % genefile
        # Get gene data from CDS file
        # Format is tab-delimited, each line has:
        # orf      chr      start     end      strand
        mapping = TabFile(genefile,column_names=('name','chr','start','end','strand'))
        n_genes_before_insert = len(gff_data)
        gff_data = GFFInsertMissingGenes(gff_data,mapping)
        print "Inserted %d missing genes" % (len(gff_data) - n_genes_before_insert)

    # Construct and insert ID for exons
    if add_exon_ids:
        print "Inserting artificial IDs for exon records"
        gff_data = GFFAddExonIDs(gff_data)

    # Construct and insert missing ID attributes
    if add_missing_ids:
        print "Inserting generated IDs for records where IDs are missing"
        gff_data = GFFAddIDAttributes(gff_data)

    # Strip attributes requested for removal
    if options.rm_attr:
        print "Removing the following attributes from all records:"
        for attr in options.rm_attr:
            print "\t* %s" % attr
        GFFUpdateAttributes(gff_data,exclude_keys=options.rm_attr)

    # Remove attributes that don't conform to KEY=VALUE format
    if strict_attributes:
        print "Removing attributes that don't conform to KEY=VALUE format"
        GFFUpdateAttributes(gff_data,exclude_nokeys=True)

    # Suppress percent encoding of attributes
    if no_attribute_encoding:
        print "Converting encoded special characters in attribute data to non-encoded form"
        logging.warning("!!! Special characters will not be correctly encoded in the output  !!!")
        logging.warning("!!! The resulting GFF may not be readable by this or other programs !!!")
        gff_data = GFFDecodeAttributes(gff_data)

    # Write to output file
    print "Writing output file %s" % outfile
    gff_data.write(outfile)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s: %(message)s")
    main()


