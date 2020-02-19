#!/bin/env python
#
#     GFTFile.py: classes for reading and manipulating data in GTF files
#     Copyright (C) University of Manchester 2012-4 Peter Briggs
#
########################################################################
#
# GFTFile.py
#
#########################################################################

"""GTFFile

Classes for reading data from GTF (Gene Transfer Format/General Transfer
Format) files.

According to the Ensembl documentation, GTF is the same as GFF version 2:

http://www.ensembl.org/info/website/upload/gff.html

however specific documentation for GTF2 can also be found at

http://mblab.wustl.edu/GTF2.html

Classes
-------

There are a number of GTF-specific classes:

 * GTFIterator: line-by-line iteration through a GTF
 * GTFFile: read data from GTF into memory so it can be easily interrogated
 * GTFDataLine: get data from a single line from a GTF file

These classes are built on top of the GFF handling classes.

Usage examples
--------------

Start by reading the data from a GTF file into a GTFFile object:

>>> gtf = GTFFile('my.gtf')

To iterate over all lines in the GTFFile object and print the
feature type:

>>> for line in gtf:
>>>    print(line['feature'])

(To iterate over all lines in the GTF without caching in memory,
printing feature types for annotation records:

>>> for line in GTFIterator('my.gtf'):
>>>    if line.type == ANNOTATION:
>>>       print(line['feature'])

)

The 'attributes' data of each annotation record are automatically converted
to a GFFAttributes object, which allows the values of named attributes to be
referenced directly. For example, if the attributes string is:

"ID=1690892742066571889;SGD=YEL026W;Gene=Sbay_5.43;Parent=Sbay_5.43"

then the value of the 'Parent' attribute can be obtained using

>>> print(line['attributes']['Parent'])

To iterate over all lines and extract 'ID' field from attributes:

>>> for line in gff:
>>>    if 'ID' in line['attributes']['ID']:
>>>       print(line['attributes']['ID'])

To iterate over all lines and print just the 'name' part of the
'ID' attribute:

>>> for line in gff:
>>>    if 'ID' in line['attributes']
>>>       print(GFFID(line['attributes']['ID']).name)

"""

from .GFFFile import GFFFile
from .GFFFile import GFFDataLine
from .GFFFile import GFFAttributes
from .GFFFile import GFFIterator
from .GFFFile import OrderedDictionary
from .GFFFile import GFF_COLUMNS

#######################################################################
# Classes
#######################################################################

class GTFDataLine(GFFDataLine):
    """Data line specific to GTF files

    Subclass of GFFDataLine specifically for handling GTF data.

    """

    def __init__(self,line=None,column_names=GFF_COLUMNS,lineno=None,
                 delimiter='\t',gff_line_type=None):
        GFFDataLine.__init__(self,line=line,column_names=column_names,
                             lineno=lineno,delimiter=delimiter,
                             gff_line_type=gff_line_type)
        self['attributes'] = GTFAttributes(str(self['attributes']))
        self._format = 'gtf'

class GTFAttributes(object):
    """Class for handling GTF 'attribute' data

    The GTF 'attribute' data consists of semi-colon separated
    data items, each of which is a 'key value' pair.

    """
    def __init__(self,attribute_data=None):
        self.__attributes = OrderedDictionary()
        self.__quotes = dict()
        for attr in GFFAttributes(attribute_data=attribute_data).nokeys():
            key = attr.split(' ')[0]
            value = ' '.join(attr.split(' ')[1:])
            if value.startswith('"') and value.endswith('"'):
                # Store quotation style for quoted values
                self.__quotes[key] = '"'
                value = value[1:-1]
            if key not in self:
                # New attribute
                self.__attributes[key] = value
            else:
                # Attribute with multiple values
                if not isinstance(self.__attributes[key],list):
                    self.__attributes[key] = [self.__attributes[key]]
                self.__attributes[key].append(value)

    def __getitem__(self,name):
        try:
            return self.__attributes[name]
        except KeyError:
            return None
    def __contains__(self,name):
        return self[name] is not None
    def __iter__(self):
        return iter(self.__attributes.keys())
    def __repr__(self):
        # Reconstruct and return the original attribute string
        attributes = list()
        for key in self:
            try:
                quotes = self.__quotes[key]
            except KeyError:
                quotes = ''
            if not isinstance(self[key],list):
                values = [self[key]]
            else:
                values = self[key]
            for value in values:
                attributes.append("%s %s%s%s;" % (key,
                                                  quotes,
                                                  value,
                                                  quotes))
        return ' '.join(attributes)

class GTFFile(GFFFile):
    """Class for handling GTF files in-memory

    Subclass of GFFFile which uses a GTFDataLine to store the
    annotation lines.

    Data from the file can then be extracted and modified using the
    methods of the GFFFile superclass (and its TabFile superclass)
    and the GTFDataLine.

    GTF is alledgedly the same as GFF version 2. See
    http://www.sanger.ac.uk/resources/software/gff/spec.html
    for the GFF specification.

    """
    def __init__(self,gtf_file,fp=None,**args):
        args['gffdataline'] = GTFDataLine
        GFFFile.__init__(self,gtf_file,fp=fp,format='gtf',**args)

class GTFIterator(GFFIterator):
    def __init__(self,gtf_file=None,fp=None,**args):
        args['gffdataline'] = GTFDataLine
        GFFIterator.__init__(self,gff_file=gtf_file,fp=fp,**args)
