#!/usr/bin/env pythons

"""
Performs an analysis of Mutual Information to the inputted fasta file.
"""

__all__ = ['mut_info']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

from os.path import expanduser
import rnatk
from rnatk.stats import Mij
from rnatk.utils import hintonDiag

def mut_info(inFile, MI_mode='free', analy_type='reduced', bioSense=False):
	"""
	Return two figures that summarize the variability and covariability in the inputted fasta file.
	Modules required:
	- rnatk
	- Mij
	- hintonDiag
	- expanduser (from os.path)
	Usage: <inFile> <MI_mode ('free' or 'basepair')> <analy_type ('reduced' or 'no_reduced')> <biological sense (default=False)>
	"""
	if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
		InFile = rnatk.fasta.openFasta(inFile)
	matrix = rnatk.fasta.transpose(inFile)
	dictResults = {}
	for indexA in range(len(matrix)):
		for indexB in range(len(matrix)):
			if bioSense:
				if indexB-indexA>3:
					if analy_type == 'reduced':
						pos1, pos2 = rnatk.fasta.basepairs(matrix[indexA], matrix[indexB])
						value = Mij(pos1, pos2, 2, mode=MI_mode)
						dictResults[indexA+1, indexB+1] = value
					else:
						value = Mij(matrix[indexA], matrix[indexB], 2, mode=MI_mode)
						dictResults[indexA+1, indexB+1] = value
			else:
				if indexB-indexA>=1:
					if analy_type == 'reduced':
						pos1, pos2 = rnatk.fasta.basepairs(matrix[indexA], matrix[indexB])
						value = Mij(pos1, pos2, 2, mode=MI_mode)
						dictResults[indexA+1, indexB+1] = value
					else:
						value = Mij(matrix[indexA], matrix[indexB], 2, mode=MI_mode)
						dictResults[indexA+1, indexB+1] = value
	hintonDiag(dictResults, expanduser('~')+'/RNA Analysis '+MI_mode+' '+analy_type+' MI_Hinton', 'RNA covariance analysis', 2)
	rnatk.fasta.info_content_file(inFile, 'RNA conservation analysis', expanduser('~')+'/RNA Analysis '+MI_mode+' '+analy_type+' Cont_Info')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("inFile", type = str, help = "input the file in FASTA format")
	parser.add_argument("--MI_mode", type = str, choices = ['free', 'basepair'], 
    	help = "How to perform the MI. free: No basepair limited to canonical. basepair: analysis limited to canonical basepairs")
	parser.add_argument("--analy_type", type = str, choices = ['reduced', 'no_reduced'],  
                        help = "reduced: only takes into account the type of basepairs presents and not the amount of these")
	parser.add_argument("--bioSense", type = bool, 
    	help = "if True is chosen, the analysis is perfomed under the restriction b-a>3, where both a and b are positions")
	args = parser.parse_args()
	haploRNA2D(args.inFile, args.MI_mode, args.analy_type, args.bioSense)
