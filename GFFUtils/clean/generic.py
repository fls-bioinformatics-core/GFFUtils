#!/usr/bin/env python
#
#     clean.generic: generic cleaning operations for GFF data
#     Copyright (C) University of Manchester 2020 Peter Briggs
#

import logging

#######################################################################
# Functions
#######################################################################

def GFFUpdateAttributes(gff_data,update_keys={},exclude_keys=[],
                        no_empty_values=True,
                        exclude_nokeys=False):
    """
    Replace and/or exclude data from the GFF attributes

    Performs manipulations on the attribute field of a GFF file,
    which typically  consists of key value pairs of the form
    'key=value', separated by semicolons.

    The operations are to update the GFF to remove (exclude)
    specific keys from the attributes of all lines, and/or
    replace the values of keys by copying from other keys in
    the same line.

    Arguments:
      update_keys: a dictionary mapping keys (attribute names)
        that should be replaced with values from other
        attributes
      exclude_keys: a list of key (attribute names) that should
        be removed from the attribute list
      no_empty_values: if set True (the default) then don't
        replace existing values with blanks (otherwise replacing
        with blank values is okay)
      exclude_nokeys: if True then any 'nokeys' attributes will
        be removed
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

def GFFAddExonIDs(gff_data):
    """
    Construct and insert a ID attribute for exons

    For each exon in the input GFF data, insert an ID
    attribute of the form:

    ID=exon:<Parent>:<n>

    where <Parent> is the name specified in the Parent
    attribute and <n> is a numerical string that is unique
    across all exons in the GFF data.

    Note that if an exon already has an ID attribute then
    it will be overwritten with the constructed ID string.

    Arguments:
      gff_data: a GFFFile object containing the GFF file
        data (which is modified in place)

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
    """
    Construct & insert a ID attribute for all features without one

    For each feature in the input GFF data that doesn't
    already have an an ID attribute, insert one of the form:

    ID=<feature>:<Parent>:<n>

    where <feature> is the feature type (e.g. 'exon', 'CDS'
    etc), <Parent> is the name specified in the Parent
    attribute and <n> is a numerical string that is unique
    across all exons in the GFF data.

    Arguments:
      gff_data: a GFFFile object containing the GFF file
        data (which is modified in place)

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
    """
    Update GFF records removing percent encoding of special characters

    The GFF format specifies that certain special
    characters in the "attribute" field must be escaped
    using "percent encoding" (i.e. HTML-style encoding).

    This function converts special characters back to
    their original form, essentially 'decoding' the
    attribute values so that they are more "human-readable".
    However the resulting attributes may be non-canonical
    and thus cause problems when fed into other GFF reading
    programs that parse the attribute field. Caveat emptor.

    Arguments:
      gff_data: a GFFFile object containing the GFF file
        data (which is modified in place)

    Returns:
      The modified GFFFile object.
    """
    for record in gff_data:
        attributes = record['attributes']
        attributes.encode(False)
    return gff_data
