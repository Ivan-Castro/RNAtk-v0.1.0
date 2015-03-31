#!/usr/bin/env python

"""
This file contains functions to analyse and build the tertiary structure
for each unique sequence in a Fasta file.

IMPORTANT!: the script 'refine_model_mmtk.py' from the module
'moderna' (version 1.7.1) must be edited. At the beginig of the file, 
in the import section add this:
from MMTK.Trajectory import StandardLogOutput

Also, change from line 233 to 242 as follow:

        energy = str(universe.energy()) #########
	print 'Initial energy = ' + energy #########
        print '\t.. starting minimization with %i cycles'%self.cycles
	minimizer = SteepestDescentMinimizer(universe, actions=[StandardLogOutput(10)]) ###############
	minimizer(steps = self.cycles)
        minimizer = ConjugateGradientMinimizer(universe, actions=[StandardLogOutput(10)]) ################
        minimizer(steps = self.cycles)
	energy = str(universe.energy()) #########
	print 'Final energy = ' + energy #########
        # write the intermediate output
        print '\t.. writing MMTK output to %s'% self.temp_pdb_file
        if self.model_passive:
            print '\t   (please note that MMTK applies a different numeration of residues.\n\t    The original one will be restored in the final output).'
        universe.writeToFile(self.mmtk_output_file)
        open(self.temp_pdb_file, 'w').write(open(self.mmtk_output_file).read())
        print '\n-------------------------------------------------------------------------------'
"""

__all__ = ['haploRNA3D']

import argparse, subprocess, time
from moderna import *
from rnatk.fasta import (openFasta, getNameSeq, seqTargTempl)
from rnatk.stats import (equLen)
from rnatk.structure3D import (replaceBase, model)
from rnatk.utils import locate

def haploRNA3D(templateInput, fastaFile, cycles, path, chain='A'):
    """
    Creates a tertiary structure for each unique sequence in a Fasta file
    using the methods from the module moderna. This function returns three
    files for each sequences: (1) a fasta file comparing the target sequence with
    the sequence of the template pdb structure, (2) a pdb file of the homologous 
    structure with no refinement, and (3) a pdb file of the homologous structure
    with refinement.
    Results are showed in a folder named 'haploRNA3D_month_day_hour_min'.
    This folder is located in Home.
    Module required:
    - argparse
    - subprocess
    - time
    - expanduser (from os.path)
    Usage: <template> <fasta file> <cycles refinement> <path> <chain (default='A')>
    """
    print 'Loading template...'
    template = load_template(templateInput, str(chain))
    clean_structure(template)
    seqTemplate = str(get_sequence(template))

    print 'Loading Fasta file...'
    records = openFasta(fastaFile)
    nameLst, seqLst = getNameSeq(records)    
    for a in nameLst:
        a = replaceBase(a, seqTemplate)
    if str(path).endswith('/'):
        nameDir = str(path)+'haploRNA3D_'+time.ctime().replace(' ','_').replace(':','_')[4:16]
    else:
        nameDir = str(path)+'/'+'haploRNA3D_'+time.ctime().replace(' ','_').replace(':','_')[4:16]
    subprocess.check_call(['mkdir', nameDir])
    muscle = locate('muscle')

    print 'Refining the models...'
    refinementScript = locate('refine_model_mmtk.py')
    secondStruct = open(nameDir+'/Secondary_Structure.fasta', 'w')
    for elm in range(len(seqLst)):
        print 'Refining', elm+1, 'of', len(seqLst)
        seqTargTempl(seqLst[elm], seqTemplate, nameLst[elm], nameDir)
        model(muscle, nameDir, nameLst[elm], templateInput, chain)
        currentTemplate = load_template(nameDir+'/'+nameLst[elm][:10])
        try:
            subprocess.check_call(['python', refinementScript, '-m',
                           nameDir+'/'+nameLst[elm][:10], '-y',
                           str(cycles), '-o', nameDir+'/'+nameLst[elm][:11]+'_refined.pdb', '-t'])
        except:
            subprocess.check_call(['python', refinementScript, '-m',
                           nameDir+'/'+nameLst[elm][:10], '-y',
                           str(cycles), '-r', '2-'+str(len(currentTemplate)),
                           '-o', nameDir+'/'+nameLst[elm][:11]+'_refined.pdb', '-t'])
        newmodel = load_model(nameDir+'/'+nameLst[elm][:11]+'_refined.pdb')
        newmodelStruc = get_secstruc(newmodel)
        newmodelSeqen = str(get_sequence(newmodel))
        secondStruct.write('>'+nameLst[elm][:11]+'\n'+newmodelSeqen+'\n'+newmodelStruc+'\n')
    secondStruct.close()
    print 'Done!'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence", type = str, help = "input the file of sequences of RNA (fasta file)")
    parser.add_argument("template", type = str, help = "input the pdb template file")
    parser.add_argument("cycles", type = str, help = "number of cycles of energy minimization")
    parser.add_argument("path", type = str, help = "Path where the folder will be created")
    parser.add_argument("--chain", type = str, help = "select the chain within the pdb template file (default = A)")
    args = parser.parse_args()
    haploRNA3D(args.template, args.sequence, args.cycles, args.path, args.chain)
