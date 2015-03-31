#!/usr/bin/env python
# coding=utf-8

"""
rnatk.fasta contains several python scripts useful to manipulate and \
analyze fasta files
"""

__all__ = ['fas2clus', 'transpose', 'openFasta', 'transcribe', 'GBcode', 'unionFasta',
           'haplotypes', 'makeFasta', 'remDupl', 'getNameSeq', 'sortConting', 'consensus',
           'purgeByKeyWord', 'positionMatch','seqTargTempl', 'noHaploSeqs', 'makeFasta_seqs_struct',
           'twoFasta', 'characterSet', 'info_content_pos', 'info_content_file', 'complementFasta',
           'basepairs', 'concernFasta', 'replaceGaps_file_byConcensus', 'fillFastaSeqs']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

import sys, re, matplotlib.pyplot as plt
from os.path import expanduser
from math import log
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import (IUPAC, Gapped)
from rnatk.seq import (seqComplem, concern, replaceGapPos_by_consensus)
from rnatk.stats import equLen

def consensus(inFile, threshold = 0.7):
    """
    Return the consensus sequences of the fasta file
    Module required:
    - AlignIO (from Bio)
    - IUPAC (from Bio.Alphabet)
    - Gapped (from Bio.Alphabet)
    - AlignInfo (from Bio.Align)
    Usage: <inFile> <threshold (default=0.7)>
    """
    alphabet = Gapped(IUPAC.ambiguous_dna)
    alignment = AlignIO.read(open(inFile), "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    return str(summary_align.dumb_consensus(threshold = threshold))

def replaceGaps_file_byConcensus(inFile, consensus):
    """
    Replace all the gaps fo each sequence by those nucleotides of the consensus that are in the same
    position.
    Module required:
    - replaceGapPos_by_consensus (from rnatk.seq)
    Usage: <inFile> <consensus>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    for index in range(len(inFile)):
        inFile[index].seq = replaceGapPos_by_consensus(str(inFile[index].seq), consensus)
    return inFile

def fillFastaSeqs(inFile, threshold=0.7, low=False, high=False, transcribe=False):
    """
    Fills the gaps of each sequence for those of the consensus sequence. Also, it enables
    to cut some section of all sequences and also for transcription.
    Module required:
    - transcribe (from rnatk.fasta)
    - concernFasta (from rnatk.fasta)
    - getNameSeq (from rnatk.fasta)
    - makeFasta (from rnatk.fasta)
    - replaceGaps_file_byConcensus (from rnatk.fasta)
    Usage: <infile> <threshold (defualt=0.7)> <low (default=False)> <high> <(default=False)> <transcribe (default=False)>
    """
    if transcribe:
        File = transcribe(inFile)
    if low and high:
        File = concernFasta(inFile, low, high)
    File = inFile
    names, seqs = getNameSeq(File)
    makeFasta(inFile, names, seqs)
    consensus_seq = consensus(File, threshold)
    File = replaceGaps_file_byConcensus(File, consensus_seq)
    names, seqs = getNameSeq(File)
    makeFasta(inFile, names, seqs)

def concernFasta(inFile, low, high):
    """
    This function does the same that the concern function but applied to the fasta file (all its sequences).
    Module required:
    - concern (from rnatk.seq)
    Usage: <inFile> <low> <high>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    for index in range(len(inFile)):
        inFile[index].seq = concern(str(inFile[index].seq), low, high)
    return inFile

def openFasta(inFile):
    """
    Return a list containing elements of Bio.SeqRecord.SeqRecord
    class.
    Modules required:
    - SeqIO (from Bio)
    Usage: <inFile>
    """
    handle = open(inFile, 'rU')
    records = list(SeqIO.parse(handle, 'fasta'))
    return records

def sortConting(inFile):
    """
    Sorts the contings from left to right.
    Modules required:
    - SeqIO (from Bio)
    Usage: <inFile>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    writing = []
    for a in range(len(inFile[0].seq)):
        memory = [] # each cycle, this list is filled with those sequence to eliminate before to run the next cycle.
        for record in inFile:
            if not str(record.seq[a]).startswith('-'):
                memory.append(record)
                writing.append(record)
        for record in memory:
            inFile.remove(record)
    return writing

def transcribe(inFile):
    """
    Transcribe the DNA to RNA for each Bio.SeqRecord.SeqRecord
    class found in the list
    Modules required:
    - SeqIO (from Bio)
    Usage: <inList>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    for a in inFile:
        a.seq = a.seq.transcribe()
    return inFile

def complementFasta(inFile, outFileName, acid='DNA', inverse=False):
    """
    Given a fasta file, this function change the sequences into its
    complement and, if inverse is True, the inverse.
    Module required:
    - seqComplem (from rnatk.seq)
    Usage: <inFile> <outFileName> <acid (default='DNA')> <inverse (default=False)>
    """
    nameLst, seqLst = getNameSeq(inFile)
    for elem in range(len(seqLst)):
        seqLst[elem] = seqComplem(seqLst[elem], acid, inverse)
    makeFasta(outFileName, nameLst, seqLst)

def fas2clus(inFile):
    """
    Convert the input file from fasta format to clustalw format.
    Modules required:
    - SeqIO (from Bio)
    Usage: <file>
    """
    SeqIO.convert(inFile, 'fasta', inFile + 'cl', 'clustal')
    clustalFile = open(inFile + 'cl')
    return clustalFile

def transpose(inFile):
    """
    This functions takes a list of sequences and compute its transpose matrix. 
    The parameter (inFile) must be a list where each element is a 
    Bio.seqRecord.seqRecord class.
    Modules required:
    - SeqIO (from Bio)
    Usage: <inFile>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    matrix = []
    for a in range(len(inFile[0].seq)):
        seq = ''
        for b in range(len(inFile)):
            seq += inFile[b].seq[a]
        matrix.append(seq)
    return matrix

def GBcode(inFile):
    """
    Change the full description of each Bio.seqRecord.seqRecord class for the
    GenBank access code only. Actually, it works for Nucleotide and SRA data base
    Modules required:
    - re
    Usage: <records> 
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    for a in inFile:
        try:
            a.description = re.findall(r'[A-Z]{1,3}[0-9]{3,9}.[0-9]{1,9}.[0-9]{1}?',a.description)[0]
        except:
            try:
                a.description = re.findall(r'[A-Z]{1,3}[0-9]{3,9}',a.description)[0]
            except:
                print 'Warning: Impossible to found a GenBank code access for %s' %(a.description)
                a.description = a.description
    return inFile

def getNameSeq(inFile):
    """
    Split the records variable into two list: nameLst
    which contains the description of each element, and
    seqLst which contains the sequence of each element.
    Module required:
    - SeqIO (from Bio)
    Usage: <inFile>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    nameLst, seqLst = [], []
    for a in inFile:
        nameLst.append(a.description)
        seqLst.append(str(a.seq))
    return nameLst, seqLst

def haplotypes(inFile): # Modify and implement the GetNameSeq function
    """
    Merge the description (or GenBank access) that share the same 
    DNA/RNA sequence. The function returns two lists: the first one
    is the description (or GenBank access), the second one is its
    respective DNA/RNA sequence.
    Module required:
    - openFasta
    Usage: <inFile>

    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    groups = {}
    a = 0
    while a < len(inFile):
        memory = []
        remove = []
        seq = str(inFile[a].seq)
        for b in inFile:
            if seq == str(b.seq):
                memory.append(b.description)
                remove.append(b)
        groups[seq]=memory
        for a in remove:
            inFile.remove(a)
        a =+ 1
    for a in groups:
        nameList = groups.values()
        seqList =  groups.keys()
    for a in range(len(nameList)):
        if len(nameList[a]) == 1:
            nameList[a] = nameList[a][0]
        else:
            name = ''
            for b in nameList[a]: # should apply here >>> name = ' '.join(nameList)
                name += ' '+b
            nameList[a] = name.strip(' ')
    return nameList, seqList

def noHaploSeqs(inFile):
    """
    This function destroys the groupment done by the haplotypes function.
    Modules required:
    - SeqIO (from Bio)
    Usage: <inFile>
    """
    saveFile = open(inFile+'_noHaploSeqs', 'w')
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    for a in range(len(inFile)):
        if ' ' in inFile[a].description:
            group = inFile[a].description.split(' ')
            for b in group:
                saveFile.write('>'+b.strip()+'\n'+str(inFile[a].seq)+'\n')
        else:
            saveFile.write('>'+inFile[a].description+'\n'+str(inFile[a].seq)+'\n')
    saveFile.close()

def makeFasta(outFile, nameList=None, seqList=None):
    """
    Make a FASTA file using two list. One list have the names, 
    the other have the sequences (or whatever).
    Usage: <nameList> <seqList> <nameFile>
    """
    FastaFile = open(str(outFile), 'w')
    if nameList and seqList:
        for a in range(len(nameList)):
            FastaFile.write('>'+nameList[a]+'\n'+seqList[a]+'\n')
    elif (bool(nameList) == False) and (bool(seqList) == True):
        for a in range(len(nameList)):
            FastaFile.write(seqList[a]+'\n')
    elif (bool(nameList) == True) and (bool(seqList) == False):
        for a in range(len(nameList)):
            FastaFile.write('>'+nameList[a]+'\n')
    else:
        print 'Warnin! No data to write to ' + str(outFile) + '. The file will be empty'
    FastaFile.close()

def makeFasta_seqs_struct(inFileSeq, inFileStruct, saveFile=False, path=expanduser('~')+'/FastaRNA_structure'):
    """
    From two fasta files, one that contains the sequences and other that contains the structures, this
    function returns a fasta-like file (dot bracket file format) which combine the sequences and the
    structures. The function have the option to save the result as a file or return the result.
    Module required:
    - openFasta
    - getNameSeq
    - expanduser (from os.path)
    Usage: <inFile_seqs> <inFile_struct> <saveFile (default=False)> <path (default='~/FastaRNA_structure')>
    """
    if (type(inFileSeq) == str) and (str(type(inFileSeq[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFileSeq = openFasta(inFileSeq)
    if (type(inFileStruct) == str) and (str(type(inFileStruct[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFileStruct = openFasta(inFileStruct)
    if len(inFileSeq) != len(inFileStruct):
        raise Exception ('Files have different sizes of sequences and structures')
    names_seqs, seqs = getNameSeq(inFileSeq)
    names_struct, struct = getNameSeq(inFileStruct)
    result = []
    for index in range(len(seqs)):
        if names_seqs[index] != names_struct[index]:
            print 'Warning! Looks like the sequence from one file does not make reference to the structure'
        result.append('>'+names_seqs[index]+'\n'+seqs[index]+'\n'+struct[index]+'\n')
    if saveFile:
        with open(path, 'w') as toSave:
            toSave.writelines(result)
    else:
        return result

def remDupl(inFile, names=False):
    """
    This function removes equal sequences leaving an only example.
    If the parameter names is True, the function only remove those
    items that are equal both in sequence and name (same identity).
    Module required:
    - SeqIO (from Bio)
    - openFasta
    - makeFasta
    Usage: <inFile> <names (default=False)>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    memory = []
    for a in range(len(inFile)):
        for b in range(len(inFile)):
            if b-a>0:
                if names:
                    if (str(inFile[a].seq) == str(inFile[b].seq)) and (inFile[a].description == inFile[b].description):
                        if inFile[b] not in memory:
                            memory.append(inFile[b])
                else:
                    if (str(inFile[a].seq) == str(inFile[b].seq)):
                        if inFile[b] not in memory:
                            memory.append(inFile[b])
    for a in memory:
        inFile.remove(a)
    return inFile

def purgeByKeyWord(inFile, keyWord, method='absent'):
    """
    Eliminate sequences that have (method = present) or not (method = absent) a specific term.
    Module required:
    - re
    - SeqIO (from Bio)
    Usage: <inFile> <keyWord> <method (default=absent)>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    memory = []
    if method == 'absent':
        for a in inFile:
            if not bool(re.search(keyWord, a.description)):
                memory.append(a)
    if method == 'present':
        for a in inFile:
            if bool(re.search(keyWord, a.description)):
                memory.append(a)
    for a in memory:
        inFile.remove(a)
    return inFile

def positionMatch(inFile, dictMatch):
    """
    Search for sequences that contains specific nucleotides at
    specific positions. These terms must be given into a dictionary.
    For example:
    dictMatch = {50:'A',51:'G',52:('C','U'),53:'U'}
    # Based on dictMatch, only is allowed sequences having the nucleotide
    A at position 50, G at position 51, C or U at position 52, and U at
    position 53.
    Module required:
    - SeqIO (from Bio)
    Scrip usage: <inFile> <dictMatch>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    memory = []
    for a in range(len(inFile)):
        for b in dictMatch:
            if not (inFile[a].seq[b-1] in dictMatch[b]):
                if not (inFile[a] in memory):
                    memory.append(inFile[a])
    for a in memory:
        inFile.remove(a)
    return inFile

def seqTargTempl(seqTarg, seqTemplate, nameSeqTarg, nameDir):
    """
    This function creates a fasta file based only in two sequences:
    >nameSeqTarg
    seqTarg
    >Template
    seqTemplate
    The fasta file is saved in nameDir
    Usage: <seqTarg> <seqTemplate> <nameSeqTarg> <path to save>
    """
    Alignment = open(nameDir+'/'+nameSeqTarg[:10]+'.fasta', 'w')
    Alignment.write('>'+nameSeqTarg[:10]+'\n'+seqTarg+'\n'+'>Template\n'+seqTemplate)
    Alignment.close()

def unionFasta(outFile, *fasta):
    """
    Merge as many fasta files as you input
    Module required:
    - openFasta
    Usage: <outFile> <fasta1, fasta2, fasta3 ... >
    """
    fasta = list(fasta)
    for index in range(len(fasta)):
        fasta[index] = openFasta(fasta[index])
    result = []
    for index in range(len(fasta)):
        for elem in fasta[index]:
            result.append('>'+elem.description+'\n')
            result.append(str(elem.seq)+'\n')
    toSave = open(outFile, 'w')
    toSave.writelines(result)
    toSave.close()

def twoFasta(inFileA, inFileB, SaveFile, action='intersection', byGB_code=False):
    """
    From 2 fasta files, the function creates a third fasta file
    based on the action:
    -intersection: Creates a third fasta file that contains only
    those sequences shared between the fasta file A and B.
    -A_diff_B: Creates a third fasta file that contains only those
    sequences present in the fasta file A but not in B.
    -symm_diff: Creates a third fasta file with those sequences that
    are member of exactly one of A and B, but no in both.
    Module required:
    - SeqIO (from Bio)
    Scrip usage: <inFileA> <inFileB> <SaveFile> <action (default='intersection')>
    """
    if (type(inFileA) == str) and (str(type(inFileA[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFileA = openFasta(inFileA)
    if (type(inFileB) == str) and (str(type(inFileB[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFileB = openFasta(inFileB)
    if byGB_code:
        inFileA = GBcode(inFileA)
        inFileB = GBcode(inFileB)
    fileResult=[]
    if action == 'intersection':
        for a in inFileA:
            for b in inFileB:
                if [a.description, str(a.seq)] == [b.description, str(b.seq)]:
                    fileResult.append('>'+a.description+'\n')
                    fileResult.append(str(a.seq)+'\n')
    if action == 'A_diff_B':
        listA = []
        listB = []
        for a in inFileA:
            listA.append([a.description, str(a.seq)])
        for b in inFileB:
            listB.append([b.description, str(b.seq)])
        for a in listA:
            if a not in listB:
                fileResult.append('>'+a[0]+'\n')
                fileResult.append(a[1]+'\n')
    if action == 'symm_diff':
        listA = []
        listB = []
        for a in inFileA:
            listA.append([a.description, str(a.seq)])
        for b in inFileB:
            listB.append([b.description, str(b.seq)])
        for a in listA:
            if a not in listB:
                fileResult.append('>'+a[0]+'\n')
                fileResult.append(a[1]+'\n')
        for a in listB:
            if a not in listA:
                fileResult.append('>'+a[0]+'\n')
                fileResult.append(a[1]+'\n')
    with open(SaveFile, 'w') as toSave:
        toSave.writelines(fileResult)

def basepairs(seq1, seq2, acid='RNA'):
    """
    Return two sequences that represents all the basepairs presents in the sequences inputted
    Usage: <sequence 1> <sequence 2>
    """
    if acid == 'RNA':
        allow = ['AA', 'AC', 'AG', 'AU',
                 'CA', 'CC', 'CG', 'CU',
                 'GA', 'GC', 'GG', 'GU',
                 'UA', 'UC', 'UG', 'UU']
    if acid == 'DNA':
        allow = ['AA', 'AC', 'AG', 'AT',
                 'CA', 'CC', 'CG', 'CT',
                 'GA', 'GC', 'GG', 'GT',
                 'TA', 'TC', 'TG', 'TT']
    total = []
    equLen(seq1, seq2)
    for index in range(len(seq1)):
        if (seq1[index]+seq2[index] in allow) and (seq1[index]+seq2[index] not in total):
            total.append(seq1[index]+seq2[index])
    pos1, pos2 = '',''
    for elem in total:
        pos1 += elem[0]
        pos2 += elem[1]
    return pos1, pos2

def characterSet(inFile):
    """
    This function return a list of all the characters present
    in the alignment or fasta file.
    Module required:
    - openFasta
    Usage: <inFile>
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    setLst = []
    for elem in inFile:
        partialSet = list(set(str(elem.seq)))
        for character in partialSet:
            if character not in setLst:
                setLst.append(character)
            else: pass
    return setLst
    
def info_content_pos(str_vector, N=None, logBase=2):
    """
    Return the information content value given one sequence
    Module required:
    - log (from math)
    Usage: <string_vector> <Number_of_characters (int)> <logBase (default = 2)>
    Reference:
    - Schneider, T. D., Stormo, G. D., & Gold, L. 1986. Information content of binding sites on
      nucleotide sequences. J. Mol. Biol. 188: 415 - 431.
    """
    NLst = list(set(str_vector))
    lenVector = float(len(str_vector))
    characFreq = {}
    obsEntropy = 0
    for charac in NLst:
        characFreq[charac] = str_vector.count(charac)/lenVector
        obsEntropy += characFreq[charac]*log(characFreq[charac],int(logBase))
    if N:
        logos = log(int(N),int(logBase)) - (obsEntropy*(-1))
    else:
        logos = log(len(N),int(logBase)) - (obsEntropy*(-1))
    return logos

def info_content_file(inFile, title, outFile, logBase=2):
    """
    Build a graphic with the logos sequences for the inputted fasta file.
    Module required:
    - matplotlib.pyplot as plt
    - log (from math)
    - transpose (from rnatk.fasta)
    Usage: <inFile> <title> <outFile> <logBase (default = 2)>
    Reference:
    - Schneider, T. D., Stormo, G. D., & Gold, L. 1986. Information content of binding sites on
      nucleotide sequences. J. Mol. Biol. 188: 415 - 431.
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    setLst = characterSet(inFile)
    transposeMatrix = transpose(inFile)
    listResult = []
    for seq in transposeMatrix:
        result = info_content_pos(seq, len(setLst), int(logBase))
        listResult.append(result)
    plt.axis(xmin=0,xmax=len(listResult),ymin=0,ymax=log((len(setLst)),int(logBase))+0.5)
    plt.title(title)
    plt.xlabel('Positions')
    plt.ylabel('Bits')
    plt.grid(True)
    plt.plot(list(range(len(listResult))), listResult, 'k.-', (0,len(listResult)),
             ((log((len(setLst)),int(logBase))),(log((len(setLst)),int(logBase)))), 'k--')
    plt.annotate('Maximum', xy=((len(listResult)/4),log((len(setLst)),int(logBase))),
                 horizontalalignment='center', xytext=((len(listResult)/4),(log((len(setLst)),int(logBase))+0.3)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=1))
    plt.text((len(listResult)/1.5), (log((len(setLst)),int(logBase))+0.2),
             '%s symbols\n%s sequences' %(len(setLst),len(inFile)))
    plt.savefig(outFile, format='svg')
    plt.clf()
