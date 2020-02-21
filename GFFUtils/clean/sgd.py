#!/usr/bin/env python
#
#     clean.sgd: cleaning operations for SGD GFF data
#     Copyright (C) University of Manchester 2020 Peter Briggs
#

import logging
from ..GFFFile import GFFID
from ..GFFFile import OrderedDictionary

#######################################################################
# Functions
#######################################################################

def GroupByID(gff_data):
    """
    Group GFF features into subsets based on ID

    Grouping is based on the ID attribute, which is assumed
    to be of the form '<prefix>:<name>:<index>' (e.g.
    'CDS:YEL0W:3').

    Consecutive features are determined to belong to the same
    gene if they have the same name and consecutive indices.

    Arguments:
      gff_data: a list of GFF data lines

    Returns:
      A list of subsets with one or more GFF data lines, where
        each subset corresponds to the same name
    """
    this_subset = []
    subsets = []
    last_index = None
    logging.debug("%d genes submitted for grouping" %
                  len(gff_data))
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

def GFFGetDuplicateSGDs(gff_data):
    """
    Return GFF data with duplicate SGD names

    Returns an OrderedDictionary where the keys are
    the duplicate SGD names and the values are
    lists of the associated GFF records, e.g.

    >>> dups = GFFGetDuplicateSGDs(gff)
    >>> for sgd in dups.keys():
    >>>    ... loops over SGD names ...
    >>>    for data in dups[sgd]:
    >>>       ... loops over duplicates ...

    Note that duplicates are determined purely by SGD
    name; no account is taken of chromosome or strand.

    Arguments:
      gff_data: a GFFFile object containing the GFF
        file data

    Returns:
      OrderedDictionary with SGDs as keys for lists
        of the duplicate TabDataLines corresponding to
        the SGD
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

def GFFResolveDuplicateSGDs(gff_data,mapping_data,duplicates,
                            overlap_margin=1000):
    """Resolve duplicate SGD names in GFF data

    Attempts to resolve duplicate SGDs by referring to a
    list of 'best' genes.

    Note that this function doesn't remove any of the data
    from the GFF input; the calling subprogram should do
    this based on the list of "discards".

    Arguments:

      gff_data: a GFFFile object containing the GFF file data
      mapping_data: a TabFile object containing candidate genes
        to insert into the GFF data if not present
      duplicates: a dictionary with keys representing SGDs
        (each key maps to a list of duplicate GFF data lines for
        that SGD) returned by the GFFGetDuplicateSGDs function
      overlap_margin: additional number of bases either side of
        candidate gene start and end positions to consider when
        looking for overlaps with duplicates

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
            logging.debug("No genes in mapping file with matching SGD to "
                          "resolve:")
            for duplicate in duplicates[sgd]:
                attr = duplicate['attributes']
                logging.debug("\t%s %s %s %s %s L%d %s" %
                              (attr['ID'],
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
                    logging.warning("Duplicate added multiple times to "
                                    "rejects list")
                rejects.append(duplicate)
        # Check if there are any matches
        if len(genes_to_duplicates.keys()) == 0:
            logging.debug("No mapping genes matched on chromosome and "
                          "strand")
            result['unresolved_sgds_no_mapping_genes_after_filter'].append(sgd)
            continue
        # Cluster duplicates for each gene and filter by overlap
        for gene in genes_to_duplicates.keys():
            # Determine overlap region
            region = (gene['start'] - overlap_margin,
                      gene['end'] + overlap_margin)
            # Group duplicates into subsets
            genes_to_duplicates[gene] = GroupByID(genes_to_duplicates[gene])
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
                            logging.warning("Duplicate added multiple "
                                            "times to rejects list")
                        rejects.append(d)
        # End of filtering process - see what we're left with
        if len(matches) == 1:
            # Resolved
            logging.debug("Duplication resolved for %s" % sgd)
            result['resolved_sgds'].append(sgd)
            for duplicate in matches[0]:
                logging.debug("\t%s %s %s %s L%d %s" %
                              (duplicate['seqname'],
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
                    logging.debug("\t%s %s %s %s L%d %s" %
                                  (duplicate['seqname'],
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
    """
    Update ID attribute of GFF data to indicate SGD groups

    For each line in the GFF data, looks for a non-blank SGD
    value in the GFF attributes and updates the ID attribute
    with the format:

    ID=CDS:<SGD>:<i>

    where <SGD> is the SGD value and <i> is an integer index
    starting from 1.

    For lines in the GFF within 5 lines of each other and
    with matching SGDs, the integer index increases by 1 each
    time to indicate that the lines form a group, for example

    CDS:YEL0W:1, CDS:YEL0W:2 etc.

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
                idx.name = sgd
                idx.index = 1
                attributes['ID'] = str(idx)
            ln = data.lineno()
            # Loop over next 5 data lines after this looking for
            # matching SGD
            for data0 in gff_data[next_ln:next_ln+5]:
                attr0 = data0['attributes']
                sgd0 = attr0['SGD']
                if sgd0 == sgd:
                    # Found a match
                    idx0 = GFFID(attr0['ID'])
                    if idx0.code != '':
                        logging.warning("ID already has code assigned "
                                        "(L%d)" % data0.lineno())
                        logging.warning("Index will be overwritten")
                    else:
                        idx0.code = "CDS"
                    idx0.name = sgd
                    idx0.index = idx.index + 1
                    attr0['ID'] = str(idx0)
                    logging.debug("%d %s\t%d %s" % (next_ln,
                                                    idx,
                                                    data0.lineno(),
                                                    idx0))
                    # Don't look any further
                    break
    # Finished grouping by SGD
    return gff_data

def GFFInsertMissingGenes(gff_data,mapping_data):
    """Insert 'missing' genes from mapping file into GFF data

    A gene is considered 'missing' from the GFF data if its
    name doesn't match any of the SGD names in the GFF data.

    Missing genes are inserted into the GFF data at the
    appropriate position based on chromosome and start position.
    
    Arguments:
      gff_data: a GFFFile object containing the GFF file data
      mapping_data: a TabFile object containing candidate genes
        to insert into the GFF data if not present
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
            logging.warning("Bad start position at L%s ('%s') "
                            "in %s, skipped" % (gene.lineno(),
                                                gene['start'],
                                                mapping_data.filename()))
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
