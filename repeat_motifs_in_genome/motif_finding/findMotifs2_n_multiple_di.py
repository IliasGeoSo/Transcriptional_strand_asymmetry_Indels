import sys
import os
import Motif_combinatorics
import mycustom 
import json
import pdb
from os.path import join
import glob
import time
from random import randint

# 640 jobs ( 1 motif per job)
jobNumber = int(sys.argv[1])
motifsPerJob = 1
lowerBound = (jobNumber-1)*motifsPerJob
upperBound = jobNumber*motifsPerJob

#The genomes to scan for the motif occurrences, finds motifs in the plus orientation of the genome
files=glob.glob("/lustre/scratch117/cellgen/team218/igs/properties/hg19/All_chr_hg19.fa")
for filed in files:
	print "This is the input argument: %s " % jobNumber
	
	# reading the sequence file
	sequenceO = mycustom.FastaFile(filed)
	sequencesL = [ i.sequence for i in sequenceO ]
	del sequenceO

	# reading the motif file
	motifO = mycustom.FastaFile("polyNs_di.fa")
	motifsL = [ i.sequence for i in motifO[lowerBound: upperBound] ]
	del motifO

	print "lowerbound is ", lowerBound
	print "upperbound is ", upperBound 

	# now finding the sequences
	result = Motif_combinatorics.findAllMotifAllSeqs(motifsL,sequencesL)

	output_dir_path = "/nfs/compgen-04/team218/ilias/indels_pancan_motif_analysis_strand_bias_polyN2/general/patient_outputs_pancan_same_1/All_chr_hg19_di"
	
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

