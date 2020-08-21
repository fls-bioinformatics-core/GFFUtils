#!/usr/bin/env python
#
# Experimental utility to remove features from a GFF file
# based on attribute values

import os
from argparse import ArgumentParser
from fnmatch import fnmatch
from GFFUtils.GFFFile import GFFIterator
from GFFUtils.GFFFile import ANNOTATION

def main():
    """
    Main program
    """
    # Process command line
    p = ArgumentParser()
    p.add_argument('gff_file',metavar="GFF_FILE",
                   help="GFF file to operate on")
    p.add_argument('-r',nargs='+',dest="remove",
                   metavar="ATTRIBUTE=PATTERN|FILE",
                   help="Features to remove; PATTERN can be an exact "
                   "match or include wildcards (e.g. 'PAC:19865399', "
                   "'PAC:1986539*'), or a file with one pattern per "
                   "line")
    p.add_argument('-o',action='store',dest='output_gff',
                   default=None,
                   help="Name of output GFF file (default is "
                   "'FILE_clean.gff')")
    args = p.parse_args()
    # Collect patterns to match to features to remove
    remove = {}
    for attr_pattern in args.remove:
        feature,pattern = attr_pattern.split('=')
        if os.path.exists(pattern):
            # 'Pattern' is a file
            # Read patterns from the file
            with open(pattern,'rt') as fp:
                for line in fp:
                    store_pattern = line.strip().strip('"')
                    try:
                        remove[feature].append(store_pattern)
                    except KeyError:
                        remove[feature] = [store_pattern]
        else:
            # Store pattern as-is
            store_pattern = pattern.strip('"')
            try:
                remove[feature].append(store_pattern)
            except KeyError:
                remove[feature] = [store_pattern]
    print(remove)
    # Output file
    if not args.output_gff:
        outbase = os.path.splitext(os.path.basename(args.gff_file))[0]
        outfile = outbase+'_clean.gff'
    else:
        outbase = os.path.splitext(os.path.basename(args.output_gff))[0]
        outfile = args.output_gff
    # Remove the features and write out
    with open(outfile,'wt') as fp:
        for line in GFFIterator(args.gff_file):
            # Flag to indicate feature should be removed
            drop_feature = False
            if line.type == ANNOTATION:
                # Get the attributes for the feature
                attrs = line['attributes']
                for attr in attrs:
                    if attr in remove:
                        # Check if the value of this attribute
                        # matches any of the stored patterns
                        # If yes then it should be dropped
                        value = attrs[attr]
                        for pattern in remove[attr]:
                            if fnmatch(value,pattern):
                                print("Dropping %s" % line)
                                print("-- %s=%s matched %s" % (attr,
                                                               value,
                                                               pattern))
                                drop_feature = True
                                break
                    if drop_feature:
                        break
            # Keep feature unless dropped
            if not drop_feature:
                fp.write("%s\n" % line)
    print("Finished")

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    main()
