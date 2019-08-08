import re,os,sys,glob
from scipy.stats import binom_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
#import seaborn as sns 
import numpy as np
import pandas as pd

revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'V':'B','H':'D','D':'H','B':'V','M':'K','K':'M','W':'W','S':'S','R':'Y','Y':'R'}[B] for B in x][::-1])

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def distance(path):
        DistancesL=[]
        DataL = reader(path)
        for v in DataL:	
		if v[-4]=="1":
			TSS = int(v[-6])+10000
			distance = -(TSS-int(v[1]))
		elif v[-4]=="-1":
			TSS = int(v[-5])-10000
			distance = TSS-int(v[1])
               	DistancesL+=[distance]
        return DistancesL

motifs =["G","GG","GGG","GGGG","GGGGG","T","TT","TTT","TTTT","TTTTT"]
colors=["pink","hotpink","orchid","mediumorchid","darkviolet","lightblue","lightskyblue","royalblue","blue","navy"]
for stk,mot in enumerate(motifs):
	# path to files with motifs, format of file is chrom start end ID strand .  chrom start end strand gene name type 
	files1  =  glob.glob(mot+"_Coordinates*plus*intersect_wo.*detail.+*upstream_10kB_downstream*")
	for dkd,i in enumerate(files1):
		motif = i.split("/")[-1].split("_")[0]
		if motif ==mot:
		
			# non-template
			plus_plus = i 
			# template
			minus_plus = i.split("Genome.bed_plus")[0]+"Genome.bed_minus"+i.split("Genome.bed_plus")[1]
			# non-template 
			minus_minus = minus_plus.split("detail.+.strand")[0]+"detail.-.strand"+minus_plus.split("detail.+.strand")[1]
			# template
			plus_minus = plus_plus.split("detail.+.strand")[0] + "detail.-.strand" + plus_plus.split("detail.+.strand")[1]
	
			DataL_plus_plus = distance(plus_plus)
			DataL_minus_minus = distance(minus_minus)
			DataL_minus_plus = distance(minus_plus)
			DataL_plus_minus = distance(plus_minus)
	
			# 
			Distances_non_template = DataL_plus_plus + DataL_minus_minus
			Distances_template = DataL_minus_plus + DataL_plus_minus
	
			bins = range(-10000,10001,100)

			hist1=np.histogram(Distances_non_template,bins=bins)[0]
			hist2=np.histogram(Distances_template,bins=bins)[0]
			Ratio=[]
			for v in range(len(hist1)):
				if float(hist2[v])!=0:
					Ratio+=[float(hist1[v])/float(hist2[v])]
					rat_bef=float(hist1[v])/float(hist2[v])
				else:
					Ratio+=[rat_bef]					
			plt.plot(range(1,len(hist1)+1,1),Ratio,label=mot,linewidth=1.2,color=colors[stk])

matplotlib.rcParams['figure.facecolor'] = 'w'
ax = plt.subplot(111)
plt.ylabel("Non-template / Template ratio of occurrences")
plt.xticks(range(1,len(bins)+1,20),[bins[st] for st in range(0,len(bins),20)])
plt.xlabel("Distance from the TSS (nt)")
plt.legend(loc="best",fontsize='small',frameon=False)
plt.ylim(0.5,1.5)
plt.xlim(0,len(bins)+1)
plt.grid()
ax.get_yaxis().tick_left()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig("test2_small_TSS.png")
plt.close()
		
