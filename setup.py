#!/usr/bin/env python
from setuptools import setup, find_packages

import re
with open('feature_merge/__version.py') as version_file:
    __versionstr__ = '.'.join(re.search(r"__version__ = \[(\d+), (\d+), (\d+)\]", version_file.read(), re.M).group()[1:])

with open('README.rst') as readme:
    setup(
        name='feature_merge',
        version=__versionstr__,
        packages=find_packages(),
        long_description=readme.read(),
        install_requires=['gffutils'],
        scripts=['bin/feature_merge'],
        url='https://github.com/brinkmanlab/feature_merge',
        license='MIT',
        author='Nolan Woods',
        author_email='nolan_w@sfu.ca',
        description='Merge overlapping features of GFF/GTF files.',
        include_package_data=True,
        test_suite='tests',
    )
