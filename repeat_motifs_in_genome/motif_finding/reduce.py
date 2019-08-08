import json
import os
from os.path import join
import glob

def loadjson(fileNameS):
	openfile = open(fileNameS)
	allD =  json.load(openfile)
	openfile.close()
	return allD

# Provide path, number of motifs (jobs) and output as json
interval = 1
numOfMotifs = 5
cellLineDir = "All_chr_hg19"
outputDir = join("patient_outputs_pancan_same_1",cellLineDir)
allD = loadjson(os.path.join(outputDir, "All_chr_hg19__outputs_"+str(interval)))
motifsL = allD.keys()

for i in range(2*interval, numOfMotifs+1, interval):
	print "currently on ", i
	tempD = loadjson(os.path.join(outputDir,"All_chr_hg19__outputs_%s"% i)) 
	allD.update(tempD)	

openfile = open("%s_same_strand_n1.json" % cellLineDir,"w")
json.dump(allD, openfile)
openfile.close()



