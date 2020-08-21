"""
Description

Setup script to install GFFUtils: utilities for working with
GFF and GTF files

Copyright (C) University of Manchester 2011-2015,2020 Peter Briggs

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
    packages = ['GFFUtils',
                'GFFUtils.clean',
                'GFFUtils.cli'],
    entry_points = { 'console_scripts': [
        'GFF3_Annotation_Extractor = GFFUtils.cli.gff3_annotation_extractor:main',
        'GFFcleaner = GFFUtils.cli.gffcleaner:main',
        'GTF_extract = GFFUtils.cli.gtf_extract:main',
        'gtf2bed = GFFUtils.cli.gtf2bed:main',]
    },
    license = 'Artistic License',
    install_requires = ['genomics-bcftbx'],
    test_suite = 'nose.collector',
    tests_require = ['nose'],
    platforms="Posix; MacOS X; Windows",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: System Administrators",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Artistic License (AFL)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    python_requires = '>=2.7',
    include_package_data=True,
    zip_safe = False
)
