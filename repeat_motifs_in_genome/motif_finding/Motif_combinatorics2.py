import re,os,sys
from collections import Counter
import numpy as np
from scipy import stats
import itertools
from sklearn import manifold # importing MDS
from sklearn.cluster import KMeans
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from random import shuffle
import copy
import pdb
from prettytable import PrettyTable
import mycustom as mc
#import pandas as pd
import operator

revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'V':'B','H':'D','D':'H','B':'V','M':'K','K':'M','W':'W','S':'S','R':'Y','Y':'R'}[B] for B in x][::-1])

def regex_mapper(s):
    """ 
        args: string input
        takes into account IUPAC code in analysis
    """
    D = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': '[AG]', 'Y': '[CT]', 'N': '[ATGC]', 'S': '[GC]', 'W': '[AT]',
         'K': '[GT]', 'M': '[AC]',
         'B': '[CGT]', 'D': '[AGT]',
         'H': '[ACT]', 'V': '[ACG]', '^': ''}
    outputString = "";
    for nucleotide in s.upper():
        outputString += D[nucleotide];
    return outputString;

def findAllMotifAllSeqs(motifL, SequencesL, include_no_matches=True):
    finalResultD = {}
    for motifS in motifL:
        finalResultD[motifS] = []
        counter = 0 
        for sequence in SequencesL:
            positionsofmatch = []
            
	    thePattern = regex_mapper(motifS)
            sequence_ = revcompl(sequence)
            thematches1 = re.finditer(thePattern, sequence_)

            for i in thematches1:
                positionsofmatch.append(len(sequence_) - i.start()-len(motifS))

            finalResultD[motifS].append(positionsofmatch)
            counter += 1
            if counter % 1000 == 0: print 'done with sequence,', counter

        print 'done with motif', motifS
    return finalResultD


def cleanMotifMatches(motif_matchesD):
    """
    args: dictionary as input
    keys: motif sequence string
    values: list of positions
    This removes Motifs that are not found in any Sequence.
    :param motif_matchesLD: This takes input of findAllMotifAllSeq. eg. {"AGAT":[], "ATTTTT":[]}
    :return:
    """
    motif_matches_no_emptyD = {}
    for motif_match in motif_matchesD.keys():

        #motif_matches_no_empty[motif_match]

        tempL = [ j[:] for j in motif_matchesD[motif_match] if j != [] ]

        if tempL != []:
            motif_matches_no_emptyD[motif_match] = copy.deepcopy(tempL)
    return motif_matches_no_emptyD


def findMotifsWithHeader(name,sequence,motifS):
	
	chrom = name
        positionsofmatch = [] #positions in consensus seq 
        positionsofmatch_ = []	#positions in revcompl of seq

        sequence=sequence.upper()
        size=len(sequence)	
        thePattern = regex_mapper(motifS)
	
        sequence_= revcompl(sequence)
        thematches1 = re.finditer(thePattern, sequence_)

        for i in thematches1:
               positionsofmatch_.append(len(sequence_) - i.start()) 

        PositionsBed=[]
        for i in positionsofmatch_:
        	if i!=[]:
        		PositionsBed+=[[chrom,int(i)-len(motifS),int(i),"RevCompl"]]
        return PositionsBed


def BedCoordinatesTopMotifs(MotifsL,SequencesL,NamesL,output_folder,window=0):

	if len(NamesL)==len(SequencesL):         #checks that the size of the headers and the sizes of the FASTA sequences is the same
	
		for i in range(len(MotifsL)):
			#if len(MotifsL[i])>3:
				BedL=[]
				for j in range(len(SequencesL)):
					BedL+=[findMotifsWithHeader(NamesL[j],SequencesL[j],MotifsL[i])]
				path_out=os.path.join(output_folder,str(MotifsL[i])+"_Coordinates_Genome.bed")
				datafile=open(path_out,"a")
				for s in BedL:
					#print s
					if s!=[]:
						for k in s:
							bed=k
							datafile.write(str(bed[0])+'\t'+str(int(bed[1])-window)+'\t'+str(int(bed[2])+window)+'\t'+str(int(bed[1])-window)+"_"+str(int(bed[2])+window)+'\t'+str(1)+'\t'+str(".")+'\n')
				datafile.close()
	return 
				

