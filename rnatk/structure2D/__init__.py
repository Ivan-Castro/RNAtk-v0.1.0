#!/usr/bin/env pythons

"""
rnatk.structure2D contains several python functions \
useful for RNA secondary structure prediction
"""

__all__ = ['Cij', 'hammDist', 'consensusRNA', 'constrStruc', 'RNAeval_free_energy',
        'secondStruct', 'secFigFile', 'fit_to_secondStruct', 'getFree_energy',
        'compare_free_energy']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

import sys, math, subprocess, re, scipy.stats, numpy
from rnatk.stats import equLen

def hammDist(str1,str2):
    """
    Find the Hamming distance between the two strings.
    Modules required:
    - sys
    - equLen
    Usage: <string 1> <string 2>
    """
    equLen(str1,str2)
    diff = 0
    for a, b in zip(str1,str2):
        if a != b:
            diff += 1
    return diff

def Cij(pos1,pos2):
    """
    Calculate the Hofacker covariation between two RNA positions
    in an alignment. Each position must be converted into a string
    of RNA.
    Modules required:
    - sys
    - hammdist
    - equLen
    Usage: <sequence1> <sequence2>
    References:
    - Hofacker et al. 2002. Secondary Structure Prediction for
      Aligned RNA Sequences. J. Mol. Biol. 319: 1059-1066.
    """
    equLen(pos1,pos2)
    allow, Cov = ['CG','GC','AU','UA','GU','UG'], 0
    for a in range(5):
        for b in range(5):
            if b>a:
                fxy, fyz = 0, 0
                for c in range(len(pos1)):
                    if pos1[c]+pos2[c] == allow[a]: fxy += 1
                    if pos1[c]+pos2[c] == allow[b]: fyz += 1
                    else: pass
                hamm = hammDist(allow[a],allow[b]) 
                Cov += (fxy/float(len(pos1)))*hamm*(fyz/float(len(pos2)))
    return Cov

def consensusRNA(inFile):
    """
    Create a RNA consensus secondary structure from a clustalw file
    of RNA aligned sequences. The function returns a list containig
    the consensus RNA sequences, the mfe structure in bracket notation,
    its energy, the free energy of the thermodynamic ensemble and the 
    frequency of the mfe structure.
    Modules required:
    -subprocess
    -RNA
    Usage: <file (clustalw format)>
    """
    output = subprocess.Popen(['RNAalifold', inFile, '-r', '--noPS'], stdout=subprocess.PIPE).stdout.readlines()
    return output

def constrStruc(seqLength, constraint):
    """
    Return a RNA secondary structure in bracket notation. The function
    requires the length of the sequence and a list or dictionary that
    contains the 2-tuple (i,j), where  i and j is represented by a '(' 
    at position i and a ')' at position j. Therefore i<j, and biologically
    j-i>3.
    Usage: <length> <list/dict-2tuple>
    """
    key = list('.'*seqLength)
    for a in constraint:
        key[a[0]-1] = '('
        key[a[1]-1] = ')'
    return ''.join(key)

def secondStruct(cons, seq):
    # TODO: When no constraint pattern is submited, make it a str of dots with equal length to the sequence.
    """
    Obtaine the secondary structure notation under the constrain
    code.
    Modules required:
    - subprocess
    Usage: <constrain pattern (dot-bracket)> <sequence>
    """
    first = subprocess.Popen(['RNAfold', '-C', '--noPS'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    second = first.communicate(input=seq+'\n'+cons)
    return (second[0][len(cons)+1:(len(cons)*2)+1])

def getFree_energy(sequence, structure=None):
    """
    This functions returns the free energy (kcal/mol) of the RNA sequence given or not some structure. The RNAfold
    used here computes automatically the best secondary structure. To get the free energy from a sequence given a strict
    restriction use the function RNAeval_free_energy from the module structure2D.
    Module required:
    - subprocess
    - RNA
    Usage: <sequence> <structure (default=None)>
    """
    first = subprocess.Popen(['RNAfold', '-C', '--noPS'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if structure == None:
        second = first.communicate(input=sequence+'\n'+str('.'*len(sequence)))
    else:
        second = first.communicate(input=sequence+'\n'+structure)
    result = re.findall('-?[0-9]+.[0-9]+', second[0])
    return float(result[0])

def RNAeval_free_energy(sequence, structure):
    """
    This functions returns the free energy (kcal/mol) of the RNA sequence given the constrained structure.
    Module required:
    - subprocess
    - RNA
    - re
    Usage: <sequence> <structure>
    """
    first = subprocess.Popen(['RNAeval'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if structure == None:
        second = first.communicate(input=sequence+'\n'+str('.'*len(sequence)))
    else:
        second = first.communicate(input=sequence+'\n'+structure)
    result = re.findall('-?[0-9]+.[0-9]+', second[0])
    return float(result[0])

def fit_to_secondStruct(secStr, sequence):
    """
    Given a RNA sequence and a RNA secondary structure in dot-bracket notation, this function try to fit the
    structure to the sequence and return such scondary structure in dot-bracket notation.
    Module required:
    - equLen (from rnatk.stats)
    Usage: <secondary structure> <RNA sequence>
    """
    equLen(secStr, sequence) # chech if the length of the structure is equal to the sequence
    charLst = list(secStr) # build a list where each item is one character of the secondary structure inputted
    if charLst.count('(') != charLst.count(')'): # this check if the amount of '(' characters is equal to the amount of
                                                 # ')' characters (which it should be!)
        print 'A bracket is missing. Quitting'
        sys.exit()
    index, search = 0, 1
    allow = ['CG','GC','AU','UA','GU','UG']
    while ('(' in charLst) and (')' in charLst):
        repeat = False
        if charLst[index] == '(': # this is the reference point (index) for search the next ')' character
            if charLst[index+search] == '.' or charLst[index+search] == 'r' or charLst[index+search] == 'l':
                search += 1
            if charLst[index+search] == '(': # if the next character is this (open-bracket), this is the new reference
                                             # point and the process restart with it.
                index += search # New reference point
                search = 1 # the search restart to 1
                repeat = True # this means that this new reference point must be ignored by the next line
            if (charLst[index+search] == ')') and repeat==False:
                if (sequence[index]+sequence[index+search]) in allow:
                    charLst[index] = 'l'
                    charLst[index+search] = 'r'
                else:
                    charLst[index] = '.'
                    charLst[index+search] = '.'
                index = 0
                search = 1
        else:
            index += 1
    result = ''.join(charLst) # here the resulted is compiled to a string variable
    return result.replace('l', '(').replace('r', ')')

def compare_free_energy(dot_bracket_file, adjust_struct):
    """
    This function compare all the RNA secondary structure's free energy obtained by the method
    used in the module rnatk.structure2D.haploRNA2D against the structure consensus propose elsewhere.
    The function returns the statistics t and the p-value plus a dictionary inside another dictionary
    with the keys "calc_struct" and "adjust_struct" each one with info about the "mean", "variance" and "mode".
    Module required:
    - RNAeval_free_energy (from rnatk.structure2D)
    - fit_to_secondStruct (from rnatk.structure2D)
    - expanduser (from os.path)
    - numpy
    - scipy.stats
    Usage: <dot bracket file> <consensus structure>
    """
    free_energy_calc_struct = []
    free_energy_adjust_struct = []
    if type(dot_bracket_file) == str:
        inFile = open(dot_bracket_file, 'r').readlines()
        for index in range(len(inFile)):
            inFile[index] = inFile[index].strip('\n').strip('>')
    if type(dot_bracket_file) == list:
        inFile = []
        for elem in dot_bracket_file:
            for item in elem.strip('\n').strip('>').split('\n'):
                inFile.append(item)
    dictResult = {}
    for index in range(len(inFile))[::3]:
        dictResult[inFile[index]] = {'seq':inFile[index+1], 'structure':inFile[index+2],
                                     'adjusted_structure':fit_to_secondStruct(adjust_struct, inFile[index+1])}
    for key in dictResult:
        free_energy_calc_struct.append(RNAeval_free_energy(dictResult[key]['seq'], dictResult[key]['structure']))
        free_energy_adjust_struct.append(RNAeval_free_energy(dictResult[key]['seq'], dictResult[key]['adjusted_structure']))
    homoVar = scipy.stats.bartlett(free_energy_calc_struct, free_energy_adjust_struct)
    if homoVar[1] < 0.05:
        mean_test = scipy.stats.ttest_ind(free_energy_calc_struct, free_energy_adjust_struct, equal_var=False)
        mean_test_performed = 'Welch\'s test'
    else:
        mean_test = scipy.stats.ttest_ind(free_energy_calc_struct, free_energy_adjust_struct, equal_var=True)
        mean_test_performed = 't-student test'
    calc_struc_mode = scipy.stats.mode(free_energy_calc_struct)
    adjust_struc_mode = scipy.stats.mode(free_energy_adjust_struct)
    calc_struc_mean = numpy.mean(free_energy_calc_struct)
    adjust_struc_mean = numpy.mean(free_energy_adjust_struct)
    calc_struc_var = numpy.var(free_energy_calc_struct)
    adjust_struc_var = numpy.var(free_energy_adjust_struct)
    return (mean_test, mean_test_performed, {'calc_struct':{'mode':calc_struc_mode, 'mean':calc_struc_mean, 'variance':calc_struc_var},
                                             'adjust_struct':{'mode':adjust_struc_mode, 'mean':adjust_struc_mean, 'variance':adjust_struc_var}})

def secFigFile(VarnaLocate, secStr, nameDir, nameFile, seq=False):
    """
    Build a secondary structure figure using Varna.
    Module required:
    - subprocess
    Usage: <path to VARNA> <secondary structure> <path to file> <nameFile> <with sequence (default=False)>
    """
    if len(nameFile)>250:
        nameFile=nameFile[0:200]+'(...)'
    if seq:
        subprocess.check_call(['java', '-cp', VarnaLocate, 'fr.orsay.lri.varna.applications.VARNAcmd', 
                          '-sequenceDBN', seq, '-structureDBN', secStr, '-o', nameDir+'/'+nameFile+'.png', 
                          '-backbone', '#000000', '-bpStyle', 'rnaviz', '-title', nameFile, 
                          '-resolution', '4.0'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        subprocess.check_call(['java', '-cp', VarnaLocate, 'fr.orsay.lri.varna.applications.VARNAcmd', 
                          '-structureDBN', secStr, '-o', nameDir+'/'+nameFile+'.png', '-backbone', '#000000', 
                          '-baseInner', '#000000', '-bpStyle', 'rnaviz', '-title', nameFile, 
                          '-resolution', '4.0'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
