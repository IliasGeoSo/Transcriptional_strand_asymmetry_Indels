import sys
import os
import Motif_combinatorics2
import mycustom 
import json
import pdb
from os.path import join
import glob
import time
from random import randint

# 32 jobs ( 1 motifs per job)
jobNumber = int(sys.argv[1])
motifsPerJob = 1
lowerBound = (jobNumber-1)*motifsPerJob
upperBound = jobNumber*motifsPerJob

#The genomes to scan for the motif occurrences, finds motifs in the minus orientation of the genome
files=glob.glob("/lustre/scratch117/cellgen/team218/igs/properties/hg19/All_chr_hg19.fa")
#read the genomic files, in this case only the human genome hg19
for filed in files:
	# reading the sequence file
	sequenceO = mycustom.FastaFile(filed)
	sequencesL = [ i.sequence for i in sequenceO ]
	del sequenceO

	# reading the motif file
	motifO = mycustom.FastaFile("polyN1.fa")
	motifsL = [ i.sequence for i in motifO[lowerBound: upperBound] ]
	del motifO

	print "lowerbound is ", lowerBound
	print "upperbound is ", upperBound 

	# now finding the sequences
	result = Motif_combinatorics2.findAllMotifAllSeqs(motifsL,sequencesL)

	#provide path to output here
	output_dir_path = "patient_outputs_pancan_opposite_1/All_chr_hg19"
	
	timer=randint(1,20)
	time.sleep(timer)
        try:
                if not os.path.exists(output_dir_path):
                        os.makedirs(output_dir_path)
        except OSError, e:
                pass

	openfile = open( join(output_dir_path,"%s_outputs_%s" % ("All_chr_hg19_", upperBound)), "w")

	json.dump(result, openfile)
	openfile.close()

