import re,os,sys,glob

datafile=open("polyNs_di.fa","w")
motifs=["AA","AT","AC","AG","TT","TA","TC","TG","CC","CG","CA","CT","GG","GT","GC","GA"]
size=range(1,11)
for s in size:
	for motif in motifs:
		datafile.write(">"+str(motif*s)+'\n')
		datafile.write(s*motif+'\n')
datafile.close()
