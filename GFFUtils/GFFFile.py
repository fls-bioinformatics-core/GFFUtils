#!/bin/env python
#
#     GFFFile.py: classes for reading and manipulating data in GFF files
#     Copyright (C) University of Manchester 2012-4 Peter Briggs
#
########################################################################
#
# GFFFile.py
#
#########################################################################

"""GFFFile

Classes for reading data from GFF (Gene-Finding Format/General Feature
Format) files (see http://www.sanger.ac.uk/resources/software/gff/spec.html
for details of the format).

Classes
-------

There are a number of GFF-specific classes:

 * GFFIterator: line-by-line iteration through a GFF
 * GFFFile: read data from GFF into memory so it can be easily interrogated
 * GFFAttributes: read data from GFF attributes field to make it easier to
   handle
 * GFFID: handle data stored in 'ID' attribute

There is an additional base class:

 * OrderedDictionary: augumented dictionary which keeps its keys in the
   order they were added to the dictionary

Usage examples
--------------

Start by reading the data from a GFF file into a GFFFile object:

>>> gff = GFFFile('my.gff')

To iterate over all lines in the GFFFile object and print the
feature type:

>>> for line in gff:
>>>    print line['feature']

(To iterate over all lines in the GFF without caching in memory,
printing feature types for annotation records:

>>> for line in GFFIterator('my.gff'):
>>>    if line.type == ANNOTATION:
>>>       print line['feature']

)

The 'attributes' data of each annotation record are automatically converted
to a GFFAttributes object, which allows the values of named attributes to be
referenced directly. For example, if the attributes string is:

"ID=1690892742066571889;SGD=YEL026W;Gene=Sbay_5.43;Parent=Sbay_5.43"

then the value of the 'Parent' attribute can be obtained using

>>> print line['attributes']['Parent']

To iterate over all lines and extract 'ID' field from attributes:

>>> for line in gff:
>>>    if 'ID' in line['attributes']['ID']:
>>>       print line['attributes']['ID']

To iterate over all lines and print just the 'name' part of the
'ID' attribute:

>>> for line in gff:
>>>    if 'ID' in line['attributes']
>>>       print GFFID(line['attributes']['ID']).name

"""

#######################################################################
# Import modules that this module depends on
#######################################################################

import logging
import copy
import urllib
from collections import Iterator
from bcftbx.TabFile import TabFile
from bcftbx.TabFile import TabDataLine

#######################################################################
# Constants/globals
#######################################################################

# Columns in GFF/GTF files
GFF_COLUMNS = ('seqname',
               'source',
               'feature',
               'start',
               'end',
               'score',
               'strand',
               'frame',
               'attributes')

# Types of line in GFF files
#
# "Pragma" line starts with '##'
PRAGMA = 0
# "Comment" line starts with single '#'
COMMENT = 1
# "Annotation" lines are tab-delimited fields containing annotation data 
ANNOTATION = 2

#######################################################################
# Class definitions
#######################################################################

class GFFDataLine(TabDataLine):
    """Data line specific to GFF files

    Subclass of TabDataLine which automatically converts data in the
    attributes column into a GFFAttributes object. This allows GFF
    attributes to be accessed directly using the syntax
    gff_line['attributes']['Parent'].
    """
    def __init__(self,line=None,column_names=GFF_COLUMNS,lineno=None,delimiter='\t',
                 gff_line_type=None):
        TabDataLine.__init__(self,line=line,column_names=column_names,
                             lineno=lineno,delimiter=delimiter)
        # Convert attributes to GFFAttributes object
        self['attributes'] = GFFAttributes(self['attributes'])
        # Metadata
        self.__type = gff_line_type

    @property
    def type(self):
        """'Type' (pragma, comment, annotation)  associated with the GFF data line

        'type' is either None, or one of the module-level constants PRAGMA, COMMENT,
        ANNOTATION, indicating the type of data held by the line.
        """
        return self.__type

class GFFFile(TabFile):
    """Class for handling GFF files in-memory

    Subclass of TabFile which uses the GFFIterator to process the
    contents of a GFF file and store annotation lines as GFFDataLines.

    Data from the file can then be extracted and modified using the
    methods of the TabFile and GFFDataLine classes.

    See http://www.sanger.ac.uk/resources/software/gff/spec.html
    for the GFF specification.
    """
    def __init__(self,gff_file,fp=None,gffdataline=GFFDataLine,format='gff'):
        # Storage for format info
        self._format = format
        self._version = None
        # Initialise empty TabFile
        TabFile.__init__(self,None,fp=None,
                         tab_data_line=GFFDataLine,
                         column_names=GFF_COLUMNS)
        # Populate by iterating over GFF file
        for line in GFFIterator(gff_file=gff_file,fp=fp,
                                gffdataline=gffdataline):
           if line.type == ANNOTATION:
                # Append to TabFile
                self.append(tabdataline=line)
           elif line.type == PRAGMA:
               # Try to extract relevant data
               pragma = str(line)[2:].split()
               if pragma[0] == 'gff-version':
                   self._version = pragma[1]

    def write(self,filen):
        """Write the GFF data to an output GFF

        Arguments:
          filen: name of file to write to
        """
        fp = open(filen,'w')
        if self.format == 'gff':
            fp.write("##gff-version 3\n")
        TabFile.write(self,fp=fp)
        fp.close()

    @property
    def format(self):
        """Return the format e.g. 'gff'

        """
        return self._format

    @property
    def version(self):
        """Return the version e.g. '3'

        """
        return self._version

class OrderedDictionary(object):
    """Augumented dictionary which keeps keys in order

    OrderedDictionary provides an augmented Python dictionary
    class which keeps the dictionary keys in the order they are
    added to the object.

    Items are added, modified and removed as with a standard
    dictionary e.g.:

    >>> d[key] = value
    >>> value = d[key]
    >>> del(d[key])

    The 'keys()' method returns the OrderedDictionary's keys in
    the correct order.
    """
    def __init__(self):
        self.__keys = []
        self.__dict = {}

    def __getitem__(self,key):
        if key not in self.__keys:
            raise KeyError
        return self.__dict[key]

    def __setitem__(self,key,value):
        if key not in self.__keys:
            self.__keys.append(key)
        self.__dict[key] = value

    def __delitem__(self,key):
        try:
            i = self.__keys.index(key)
            del(self.__keys[i])
            del(self.__dict[key])
        except ValueError:
            raise KeyError

    def __len__(self):
        return len(self.__keys)

    def __contains__(self,key):
        return key in self.__keys

    def __iter__(self):
        return iter(self.__keys)

    def keys(self):
        return copy.copy(self.__keys)

    def insert(self,i,key,value):
        if key not in self.__keys:
            self.__keys.insert(i,key)
            self.__dict[key] = value
        else:
            raise KeyError, "Key '%s' already exists" % key

class GFFAttributes(OrderedDictionary):
    """Class for handling GFF 'attribute' data

    The GFF 'attribute' data consists of semi-colon separated
    data items, which might be single values or 'key=value'
    pairs.

    Given an input string of this form, a GFFAttributes
    object provides access to the keys via the keys() method
    (lists all keys in the order they were encountered), and
    access to values via the GFFAttributes[key] syntax.

    Keyed data can be modified and deleted via this syntax.

    The values without keys can be accessed as a list via the
    nokeys() method. Special characters that have been escaped
    using "percent encoding" (e.g. semicolons, equals, tabs etc)
    are automatically converted back to their original values
    (e.g. ';' rather than '%3B').

    str(GFFAttributes) returns the attribute data as a
    string representation which restores the original format
    found in the GFF file (i.e. semi-column separated data
    items with single values or 'key=value' pairs as
    appropriate). Values which contain special characters
    will be escaped appropriately using URL percent encoding
    (i.e. the reverse of the decoding process).
    """
    def __init__(self,attribute_data=None):
        OrderedDictionary.__init__(self)
        self.__nokeys = []
        # Special attributes which can have multiple values
        self.__multivalued_attributes = ('Parent',
                                         'Alias',
                                         'Note',
                                         'Dbxref',
                                         'Ontology_term')
        # Flag indicating whether to encode values on output
        self.__encode_values = True
        # Flag indicating whether data came with trailing semicolon
        self.__trailing_semicolon = False
        # Extract individual data items
        if attribute_data:
            for item in attribute_data.split(';'):
                if not item:
                    continue
                try:
                    i = item.index('=')
                    key = item[:i].strip()
                    value = item[i+1:].strip()
                except ValueError:
                    # No equals sign
                    key = ''
                    value = item.strip()
                # Percent-decode value
                value = urllib.unquote(value)
                # Store data
                if key == '':
                    # No key: store in a list
                    self.__nokeys.append(value)
                else:
                    # Store key-value pair
                    self[key] = value
            self.__trailing_semicolon = attribute_data.endswith(';')

    def nokeys(self):
        return self.__nokeys

    def encode(self,new_setting=None):
        """Query or set the value of the 'encoding' flag

        By default attribute values are encoded when written as a
        string i.e. special characters are escaped using HTML-style
        encoding to convert them to "percent codes" e.g. '%3B'.
        This method allows this encoding to be turned on or off
        explicitly.

        To query the current setting, invoke the method without
        arguments e.g.:

        >>> attr.encode()
        ... True

        To turn encoding on or off, invoke with a boolean argument,
        e.g.:

        >>> attr.encode(False)  # Turns encoding off
        ... False

        Note that turning the decoding off is not recommended when
        writing the attributes back to a GFF file.
        
        Arguments:
          new_setting: (optional) the new value for the encoding
            flag - True to turn encoding on, False to turn it off

        Returns:
          Value of the encoding flag, indicating whether encoding
          is on (True) or off (False).
        """
        if new_setting in (True,False):
            self.__encode_values = new_setting
        elif new_setting is not None:
            raise ValueError, "bad value '%s': can only be True or False" % new_setting
        return self.__encode_values

    def __escape_value(self,key,value):
        """Internal: return escaped value of input string

        Performs the URL percent encoding of values for the GFF
        attributes.

        Note that for certain predefined keys, an unescaped comma
        is allowed in the value. For all others commas will be
        escaped.

        Arguments:
          key: name of the attribute that the value belongs to
          value: the string to be encoded
        """
        if key in self.__multivalued_attributes:
            return urllib.quote(value,safe=" ,:^*$@!+?|")
        else:
            return urllib.quote(self[key],safe=" :^*$@!+?|")

    def __repr__(self):
        items = []
        if self.__encode_values:
            # Percent encode the attributes
            for key in self.keys():
                items.append("%s=%s" % (key,self.__escape_value(key,self[key])))
        else:
            # Don't encode the attributes
            for key in self.keys():
                items.append("%s=%s" % (key,self[key]))
        # Add nokeys items
        for item in self.__nokeys:
            items.append(item)
        # Dealing with trailing semicolon
        if self.__trailing_semicolon:
            items.append('')
        return ';'.join(items)

class GFFID(object):
    """Class for handling ID attribute in GFF data

    The GFFID class takes an ID string which is expected to be of
    one of the forms

    <name> (e.g. 'XYZ123-A')

    or
    
    <code>:<name>:<index> (e.g. 'CDS:XYZ123-A:2')

    The components can be queried and changed via the code,
    name and index properties of the class, e.g.

    >>> code = GFFID('CDS:XYZ123-A:2').code

    etc

    If there is no code then this is returned as an empty string;
    if there is no index then this is returned as 0.

    Doing e.g. str(gffid) returns the appropriate string
    representation.
    """
    def __init__(self,gff_id):
        items = str(gff_id).split(':')
        if len(items) == 1:
            self.code = ''
            self.name = gff_id
        else:
            self.code = items[0]
            self.name = items[1]
        if len(items) > 2:
            try:
                self.index = int(items[2])
            except ValueError:
                print "Bad id? %s %s" % (gff_id,items)
                raise ValueError
        else:
            self.index = 0
            
    def __repr__(self):
        if not self.code:
            return "%s" % self.name
        else:
            return "%s:%s:%d" % (self.code,self.name,self.index)

class GFFIterator(Iterator):
    """GFFIterator

    Class to loop over all records in a GFF file, returning a GFFDataLine
    object for each record.

    Example looping over all reads
    >>> for record in GFFIterator(gff_file):
    >>>    print record
    """

    def __init__(self,gff_file=None,fp=None,gffdataline=GFFDataLine):
        """Create a new GFFIterator

        Arguments:
           gff_file: name of the GFF file to iterate through
           fp: file-like object to read GFF data from
           gffdataline: GFFDataLine-like class to instantiate
             and return for each record in the GFF
        """
        if fp is not None:
            self.__fp = fp
            self.__close_fp = False
        else:
            self.__fp = open(gff_file,'rU')
            self.__close_fp = True
        self.__gffdataline = gffdataline
        self.__lineno = 0

    def next(self):
        """Return next record from GFF file as a GFFDataLine object
        """
        line = self.__fp.readline()
        self.__lineno += 1
        if line != '':
            # Set type for line
            if line.startswith("##"):
                # Pragma
                type_ = PRAGMA
            elif line.startswith("#"):
                # Comment line
                type_ = COMMENT
            else:
                # Annotation line
                type_ = ANNOTATION
            # Convert to GFFDataLine
            return self.__gffdataline(line=line,lineno=self.__lineno,gff_line_type=type_)
        else:
            # Reached EOF
            if self.__close_fp: self.__fp.close()
            raise StopIteration
