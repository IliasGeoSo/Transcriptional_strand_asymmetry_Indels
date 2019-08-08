import re,os,sys,glob,json
from collections import Counter
import Motif_combinatorics
import Motif_combinatorics2
import numpy as np
from scipy import stats
import itertools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
from random import shuffle
import copy
from prettytable import PrettyTable
import mycustom as mc
#import pandas as pd
import operator
from scipy.stats import binom_test

revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'V':'B','H':'D','D':'H','B':'V','M':'K','K':'M','W':'W','S':'S','R':'Y','Y':'R'}[B] for B in x][::-1])

if 1==1:
	# provide path to the two json files generated for each strand
	json_file1 = "All_chr_hg19_di_same_strand_n1.json"
	json_file2 = "All_chr_hg19_di_opposite_strand_n1.json"

	# provide path to genome
	seqFileName1 = "/lustre/scratch117/cellgen/team218/igs/properties/hg19/All_chr_hg19.fa"
	seqFileName2 = "/lustre/scratch117/cellgen/team218/igs/properties/hg19/All_chr_hg19.fa"
    
	fastaO = mc.FastaFile(seqFileName1)
	SequencesL1 = [i.sequence for i in fastaO.getSequences()]
	NamesL1 = [i for i in fastaO.getNames()]
	sequence_length1 = len(SequencesL1[0]) 
	nucs_seq1 = sequence_length1*len(SequencesL1)

	fastaO = mc.FastaFile(seqFileName2)
	SequencesL2 = [i.sequence for i in fastaO.getSequences()]
	NamesL2 = [i for i in fastaO.getNames()]
	sequence_length2 = len(SequencesL2[1])
	nucs_seq2 = sequence_length2*len(SequencesL2)

	# provide file of polyN motifs (e.g. di-nucleotides or mono-nucleotides)
	motifO = mc.FastaFile("polyNs_di.fa")
	consensusL = [i.sequence for i in motifO.getSequences()]
	MotifsL = consensusL

	if 'motif_matchesD' not in globals().keys():
		openfile = open(json_file1)
	        motif_matchesD1 = json.load(openfile)
        	openfile.close()
	
	if 'motif_matchesD' not in globals().keys():
		openfile = open(json_file2)
		motif_matchesD2 = json.load(openfile)
		openfile.close()


	MotifsL2=[]
	for i in MotifsL:
		if i not in MotifsL2 and revcompl(i) not in MotifsL2:
			MotifsL2+=[i]

        Motifs_Enriched=[]
        for motif in MotifsL:
                if motif not in Motifs_Enriched:
                        Motifs_Enriched+=[motif]
	
	# generate folders to put the output files
	os.system("mkdir "+ "All_chr_hg19_TSA_di")
	os.system("mkdir "+ "All_chr_hg19_TSA_di/plus/")
	os.system("mkdir "+ "All_chr_hg19_TSA_di/minus/")
	Motif_combinatorics.BedCoordinatesTopMotifs(Motifs_Enriched,SequencesL1,NamesL1,"All_chr_hg19_TSA_di/plus/",window=0)
	Motif_combinatorics2.BedCoordinatesTopMotifs(Motifs_Enriched,SequencesL2,NamesL2,"All_chr_hg19_TSA_di/minus/",window=0)


