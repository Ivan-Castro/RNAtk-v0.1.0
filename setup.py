#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

readme = open('README.txt').read()

setup(
    name='rnatk',
    version='0.1.0',
    author='Ivan Castro',
    author_email='sivanc7@gmail.com',
    maintainer='Ivan Castro',
    description='RNA toolkit',
    long_description=readme,
    classifiers=[
         'Development Status :: 1 - Planning',
         'Environment :: Console',
         'Intended Audience :: Science/Research',
         'License :: OSI Approved :: GNU General Public License (GPL)',
         'Natural Language :: English',
         'Operating System :: OS independent',
         'Programming Language :: Python :: 2.7',
         'Topic :: Scientific/Engineering :: Bio-informatics',
         'Topic :: Software Development :: Libraries :: Python Modules',
         ],
    platforms=['Linux'],
    license=['GPL'],
    keywords=['RNA', 'Bioinformatics', 'Biology', 'Structure'],
    packages=find_packages(),
    requires=['RNA (>=2.0)', 'Bio (>=1.63)', 'numpy (>=1.8.2)', 'scipy (>=0.13.3)',
              'matplotlib (>=1.3.1)', 'moderna (>=1.7.1)', 'ete2 (>=2.2)', 'pylab',
              'networkx (>=1.9.1)']
    )
