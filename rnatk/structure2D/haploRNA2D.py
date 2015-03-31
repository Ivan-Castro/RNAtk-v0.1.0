#!/usr/bin/env python

"""
This file contains functions to analyse the secondary structure
for each unique sequence in a Fasta file.
"""

__all__ = ['haploRNA2D']

import argparse, time, subprocess
from os.path import expanduser
from rnatk.fasta import (openFasta, fas2clus, GBcode, haplotypes, transpose, makeFasta, transcribe, info_content_file)
from rnatk.seq import randomSeq
from rnatk.stats import (kernelXij, Mij, searchXijc, fig_plot, equLen_file)
from rnatk.structure2D import (Cij, consensusRNA, constrStruc, secondStruct, secFigFile)
from rnatk.utils import (locate, retrieve_varName,  
     bestPartner, bestPartner, signXij, hintonDiag, dictXij)
from rnatk.dendro import (dendroNJ, saveTree)

def haploRNA2D(inFile, title, repeats, siglevel, bootstrap, trc=False):
    """
    Creates a secondary structure for each unique sequence and a
    consensus too. All the process involves the MFE along covariance
    analysis. It also creates Fasta-like files summarizing the
    results.
    Results are showed in a folder named 'haploRNA2D_month_day_hour_min'.
    This folder will be available in Home.
    Module required:
    - time
    - expanduser (from os.path)
    Usage: <Fasta file> <title> <repeats> <significance level> <bootstrap> <transcription (default=False)>
    """
    print 'Reading file...'
    records = GBcode(inFile) # this step change the file to Bio.SeqIO format and also take only the GenBank code as reference
    equLen_file(records)
    if trc == True:
        print 'Transcribing sequences...'
        records = transcribe(records)
    fas2clus(inFile) # make a clustal file to be used by the consensus function
    consensusLst = consensusRNA(inFile+'cl')
    if consensusLst:
        consensusFile = True
        lenSeq = len(consensusLst[0])-1
        consensuStr = consensusLst[1][0:lenSeq]
    else:
        consensusFile = False
        print 'Warning! No consensus structure defined!'
        lenSeq = len(str(records[0].seq))
    nameDir = expanduser('~')+'/'+'haploRNA2D_'+time.ctime().replace(' ','_').replace(':','_')[4:16]

    print 'Making folder...'
    subprocess.check_call(['mkdir', nameDir]) # make a folder inside the home folder as nameDir shows

    VarnaLocate = locate('VARNA') # search and save the path where the application VARNA is located into the variable 'VarnaLocate' 
    if consensusFile:
        print 'Building RNA secondary structure...'
        ConsensusFasta = open(nameDir+'/'+'Consensus_Structure', 'w')
        ConsensusFasta.write('>Consensus secondary structure\n'+consensuStr)
        ConsensusFasta.close()
        secFigFile(VarnaLocate, consensuStr, nameDir, 'Consensus_Secondary_Structure')

    print 'Calculating random covariances...'
    MijLst, CijLst = [], []
    for a in range(repeats):
        seq2Tpl = (randomSeq(len(records), 'RNA'), randomSeq(len(records), 'RNA')) # make a 2-tuple of random sequences
        MijLst.append(Mij(seq2Tpl[0],seq2Tpl[1],2)) # compute and store the Mij value of the two random sequences
        CijLst.append(Cij(seq2Tpl[0],seq2Tpl[1])) # compute and store the Cij value of the two random sequences

    print 'Calculating real covariances...'
    Mijc = searchXijc(kernelXij(MijLst), siglevel, -0.5, 2.0)[0] # here is defined the Mij critical value
    Cijc = searchXijc(kernelXij(CijLst), siglevel, -0.5, 0.75)[0] # here is defined the Cij critical value
    transposeRecords = transpose(records) # here, the alignment matrix is transposed and...
    dicMij, dicCij = dictXij(transposeRecords, 'Mij'), dictXij(transposeRecords, 'Cij') # ... its Mij and Cij values are calculated here

    print 'Saving Hinton diagrams...'
    hintonDiag(dicMij, nameDir+'/Mij Hintong', 'Hinton diagram for Mij values', 2)
    hintonDiag(dicMij, nameDir+'/Cij Hintong', 'Hinton diagram for Cij values', 1.5)

    print 'Saving png KDE images for each array...'
    dicMijVal = dicMij.values()
    dicCijVal = dicCij.values()
    variables = [MijLst, CijLst, dicMijVal, dicCijVal]
    for a in range(len(variables)):
        varName = retrieve_varName(variables[a])[0]
        fig_plot(variables[a], -1, 3, nameDir+'/'+varName, varName)

    print 'Saving svg Logos images for the alignment...'
    info_content_file(records, title, nameDir+'/Logos')

    print 'Searching significative values and best partner...'
    dicMij, dicCij = signXij(dicMij, Mijc), signXij(dicCij, Cijc) # here are saved only those Mij and Cij whose value is greater than the critical value
    dicMijSign, dicCijSign = dicMij.copy(), dicCij.copy()
    dicMij, dicCij = bestPartner(dicMij), bestPartner(dicCij) # only those pair positions with higgest values are saved

    print 'Building the constraint pattern...'
    constraint = []
    for a in dicMij: # 'a' is a 2-tuple which point what positions are basepaired
        if a in dicCij: # this means: if 'a' is in both dicMij and dicCij...
            constraint.append(a) # ... store 'a' in constraint

    print 'Creating files...'
    nameFasta, seqFasta = haplotypes(records)
    makeFasta(nameDir+'/Haplotypes_Fasta', nameFasta, seqFasta)
    conStruct = constrStruc(lenSeq, constraint) # makes the secondary structure pattern for constraint
    haploStruct = []
    for a in seqFasta:
        resultHaploStruct = secondStruct(conStruct, a) 
        haploStruct.append(resultHaploStruct)
    makeFasta(nameDir+'/Haplotype_Structures_Fasta', nameFasta, haploStruct)
    others = open(nameDir+'/Other_Results', 'w')
    others.write('Critical value for Mij:\n'+str(Mijc)+'\n'+
                 'Critical value for Cij:\n'+str(Cijc)+'\n'+
                 'Significative values for Mij:\n'+str(dicMijSign)+'\n'+
                 'Significative values for Cij:\n'+str(dicCijSign)+'\n'+
                 'Mij values after "bestPartner":\n'+str(dicMij)+'\n'+
                 'Cij values after "bestPartner":\n'+str(dicCij)+'\n'+
                 'Final constraint basepairs:\n'+str(constraint)+'\n'+
                 'Consensus constraint:\n'+conStruct+'\n')
    if not consensusFile:
        others.write('No consensus structure defined!')
    others.close()

    print 'Saving phylogenetic tree of secondary structures...'
    treeNJ = dendroNJ(nameDir+'/Haplotype_Structures_Fasta', replicate=bootstrap)
    saveNJMatrix = open(nameDir+'/NJMatrix','w')
    saveNJMatrix.write(treeNJ)
    saveNJMatrix.close()
    saveTree(treeNJ, nameDir+'/Phylo_Secondary_Structure', ladderize=True)

    print 'Making png secondary structures for each haplo-sequence...'
    for a in range(len(nameFasta)):
        print 'Making',a+1,'of', len(nameFasta)
        if '(' in haploStruct[a]:
            secFigFile(VarnaLocate, haploStruct[a], nameDir, nameFasta[a], seqFasta[a])
        else:
            print nameFasta[a], 'has no a RNA secondary structure'

    print 'Done!'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type = str, help = "input the file in FASTA format")
    parser.add_argument("title", type = str, help = "Title for Logos image")
    parser.add_argument("repeats", type = int, help = "Number of repetitions used for calculate the PDF of Mij and Cij")
    parser.add_argument("siglevel", type = float, choices = [0.01, 0.05, 0.10, 0.15],  
                        help = "Significance level to determinate the covariation between pair of positions")
    parser.add_argument("bootstrap", type = int, help = "Bootstrap for the phylogenetic tree of secondary structures")
    parser.add_argument("--trc", action = 'store_true', help = "Transcribe the DNA sequence to RNA")
    args = parser.parse_args()
    haploRNA2D(args.input, args.title, args.repeats, args.siglevel, args.bootstrap, args.trc)
