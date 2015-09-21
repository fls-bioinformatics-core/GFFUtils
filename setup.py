"""
Description

Setup script to install GFFUtils: utilities for working with
GFF and GTF files

Copyright (C) University of Manchester 2011-2015 Peter Briggs

"""

readme = open('README.rst').read()

# Setup for installation etc
from setuptools import setup
import GFFUtils
setup(
    name = "GFFUtils",
    version = GFFUtils.get_version(),
    description = "utilities for working with GFF and GTF files",
    long_description = readme,
    url = 'https://github.com/fls-bioinformatics-core/GFFUtils',
    maintainer = 'Peter Briggs',
    maintainer_email = 'peter.briggs@manchester.ac.uk',
    packages = ['GFFUtils',],
    entry_points = { 'console_scripts': [
        'GFF3_Annotation_Extractor = GFFUtils.GFF3_Annotation_Extractor:main',
        'GFFcleaner = GFFUtils.GFFcleaner:main',
        'GTF_extract = GFFUtils.GTF_extract:main',]
    },
    license = 'Artistic License',
    install_requires = ['genomics'],
    test_suite = 'nose.collector',
    tests_require = ['nose'],
    include_package_data=True,
    zip_safe = False
)
