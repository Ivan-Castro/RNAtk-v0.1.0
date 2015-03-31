#!/usr/bin/env python

"""
rnatk.seq contains several python fucntions useful to DNA/RNA sequences \
manipulation
"""

__all__ = ['randomSeq', 'concern', 'seqComplem', 'replaceGapPos_by_consensus']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

import sys, random, string
from rnatk.stats import equLen

def randomSeq(lenghtSeq, acid='RNA'):
    """
    Creates a random DNA/RNA sequences of length=lenghtSeq.
    Each nucleotide at some position has equal probability
    (25%).
    Module required:
    - random
    Usage: <lenghtSeq> <acid=('DNA' or 'RNA')>
    """
    seq = ''
    for b in range(int(lenghtSeq)):
        num = random.randint(1,4)
        if num == 1: seq += 'A'
        if num == 2: seq += 'C'
        if num == 3:
            if acid=='DNA': seq += 'T'
            else: seq += 'U'
        if num == 4: seq += 'G'
    return seq

def concern(seq, low, high):
    """
    Take a section of the sequence from 'low' until 'high'.
    Usage: <sequence> <from> <until>
    """
    seq = seq[int(low):int(high)]
    return seq

def replaceGapPos_by_consensus(sequence, consensus):
    """
    Give a sequence, this function replace the gaps by the correspondent
    nucleotide from the consensus sequence.
    Module required:
    - equLen (rnatk.stats)
    Usage: <sequence> <consensus>
    """
    equLen(sequence, consensus)
    sequence = list(sequence)
    for index in range(len(sequence)):
        if sequence[index] == '-':
            sequence[index] = consensus[index]
    return ''.join(sequence)

def seqComplem(strand, acid='DNA', inverse=False):
    """
    Return the complementary sequence of a DNA/RNA sequence and
    optionally the inverse orientation
    Module required:
    - string
    Usage: <sequence> <acid=('DNA' or 'RNA')> <inverse (dafault=False)
    """
    if acid == 'DNA':
        table = string.maketrans('TAGCtagc', 'ATCGATCG')
    if acid == 'RNA':
        table = string.maketrans('UAGCuagc', 'AUCGAUCG')
    sequence = strand.translate(table)
    if inverse == True:
        sequence = sequence[::-1]
    return sequence
