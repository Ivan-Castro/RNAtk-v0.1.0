#!/usr/bin/env python

"""
rnatk is a module focused to manipulate sequences, fasta files \
and RNA secondary - tertiary structure
"""

import sys

__all__ = ['dendro', 'fasta', 'seq', 'stats', 'structure2D', 'structure3D', 
           'utils']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, rnatk Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

#Check python version
if sys.version_info < (2, 7):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError("Python-2.7 or greater is required, Python-%s used." % py_version)
