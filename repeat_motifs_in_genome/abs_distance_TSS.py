import re,os,sys,glob
from scipy.stats import binom_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
import seaborn as sns 
import numpy as np

# path to polyN files e.g.  polyTs and polyGs for plus and minus orientation, replace with the path to your motif files in plus and minus genom
files=glob.glob("../plus/T_*bed")+glob.glob("../minus/T_*bed")
# path to gene files separated by plus and minus orientation with 10kB upstream and downstream further from the TSS and TES respectively. 
files2=["genes.detail.-.strand_10kB_upstream_10kB_downstream.txt","genes.detail.+.strand_10kB_upstream_10kB_downstream.txt"]
#sort files to perform bedtools command
for v in files:
	os.system("sort -k1,1 -k2,2n " +v + " > " + v +".sorted")
	os.system("mv "+v +".sorted " +v)

for v in files:
	for j in files2:
		if "/plus/" in v:
        		os.system("bedtools intersect -a " + v + " -b " + j + " -wo -sorted > " + v.split("/")[-1]+"_plus.intersect_wo."+j.split("/")[-1])
		else:
			os.system("bedtools intersect -a " + v + " -b " + j + " -wo -sorted > " + v.split("/")[-1]+"_minus.intersect_wo."+j.split("/")[-1])


