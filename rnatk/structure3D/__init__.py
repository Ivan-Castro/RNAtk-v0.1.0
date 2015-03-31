#!/usr/bin/env pythons

"""
rnatk.structure3D contains several python functions \
useful for RNA tertiary structure prediction
"""

__all__ = ['replaceBase', 'model']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

import sys, subprocess
from moderna import *
from rnatk.fasta import (openFasta, getNameSeq)
from rnatk.stats import (equLen)

def replaceBase(seqTarg, seqTemplate):
    """
    This file replaces the ambigous nucleotide 'N' from 
    a determinate position from the target sequence for
    the nucleotide X from the sequence template that belongs
    to the same position.
    Usage: <sequence Target> <sequence template>
    """
    lstTarg = list(seqTarg)
    returnStr = ''
    for a in range(len(lstTarg)):
        if lstTarg[a] == 'N':
            lstTarg[a] = seqTemplate[a]
        returnStr += lstTarg[a]
    return returnStr

def model(muscle, nameDir, nameSeqTarg, template, chain='A'):
    """
    Creates a pdb file with a structure homologous to the
    template. It needs a fasta file that contains only the
    (1) target sequence and (2) the template sequence.
    Module required:
    - subprocess
    - * (from moderna)
    Usage: <path to MUSCLE> <path to File> <name file> <pdb template> <chain (default='A')>
    """
    try:
        subprocess.check_call([muscle, '-in', nameDir+'/'+nameSeqTarg[:10], 
            '-out', nameDir+'/'+nameSeqTarg[:10]])
        alignModel = load_alignment(nameDir+'/'+nameSeqTarg[:10])
    except:
        subprocess.check_call([muscle, '-in', nameDir+'/'+nameSeqTarg[:10]+'.fasta', 
            '-out', nameDir+'/'+nameSeqTarg[:10]+'.fasta'])
        alignModel = load_alignment(nameDir+'/'+nameSeqTarg[:10]+'.fasta')
    templ = load_template(template, chain)
    model = create_model(templ, alignModel)
    model.write_pdb_file(nameDir+'/'+nameSeqTarg[:10])
