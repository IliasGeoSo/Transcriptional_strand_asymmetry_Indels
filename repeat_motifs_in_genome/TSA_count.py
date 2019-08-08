import re,os,sys,glob
from scipy.stats import binom_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# provide path of polyNs
files=glob.glob("All_chr_hg19_TSA/plus/*bed")#+glob.glob("All_chr_hg19_TSA/minus/*bed")
genes_plus = 'genes.detail.+.strand.txt'
genes_minus = 'genes.detail.-.strand.txt'
for v in files:
	os.system("bedtools intersect -a " + v + " -b " + genes_plus + " -u > " + v+".genes_plus_strand")
	os.system("bedtools intersect -a " + v + " -b " + genes_minus + " -u > " + v + ".genes_minus_strand")


