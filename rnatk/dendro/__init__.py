#!/usr/bin/env python

"""
rnatk.dendro contains several python function useful to manipulate \
phylogenetics data
"""

__all__ = ['distance', 'dendroNJ', 'saveTree']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

import sys
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
from Bio.Phylo.Consensus import *
from ete2 import Tree, TreeStyle

def distance(inFile, model='identity'):
    """
    Given an alingment file (in fasta format), this function return a distance matrix.
    Module required:
    - AlignIO (from Bio)
    - DistanceCalculator (from Bio.Phylo.TreeConstruction)
    Usage: <inFile> <model (default = 'identity')>
    """
    aln = AlignIO.read(inFile, 'fasta') # read the alignment
    calculator = DistanceCalculator(model) # prepare the mode to calculate the distance
    dm = calculator.get_distance(aln) # calculate the distance of the alignment
    return dm

def dendroNJ(inFile, model='identity', bootstrap=True, replicate=100):
    """
    Given an alingment in fasta format, the function returns a Neighbor Joining tree in newick format.
    Module required:
    - AlignIO (from Bio)
    - DistanceCalculator (from Bio.Phylo.TreeConstruction)
    - DistanceTreeConstructor (from Bio.Phylo.TreeConstruction)
    - bootstrap_consensus (from Bio.Phylo.Consensus)
    Usage: <inFile> <model (default = 'identity')> <bootstrap (default = True)>
                           <replicate (default = 100)>
    """
    aln = AlignIO.read(inFile, 'fasta') # read the alignment
    constructor = DistanceTreeConstructor(DistanceCalculator(model), 'nj')
    if bootstrap:
        tree = bootstrap_consensus(aln, int(replicate), constructor, majority_consensus)
    else:
        tree = constructor.build_tree(aln)
    return tree.format('newick')

def saveTree(inFile, fileSave, formatFile='svg', ladderize=False):
    """
    Given a file containing a tree in newick format, the function
    save an imagen of such tree.
    Module required:
    - Tree (from ete2)
    - TreeStyle (from ete2)
    Usage: <inFile> <path/name to save> <formatFile (default='svg')>
                           <ladderize (default=False)>
    """
    tree = Tree(inFile)
    if ladderize:
        tree.ladderize()
    ts = TreeStyle()
    ts.show_branch_support = True
    tree.render(fileSave+'.'+formatFile, tree_style=ts)
