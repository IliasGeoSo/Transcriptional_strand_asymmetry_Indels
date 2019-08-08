import re,os,sys,glob
from scipy.stats import binom_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
import seaborn as sns 
import numpy as np
from random import randint
import time

files1=glob.glob("All_chr_hg19_TSA/plus/*bed")
files2=glob.glob("All_chr_hg19_TSA/minus/*bed")
for v in files1:
                print v
                os.system("sort -k1,1 -k2,2n " +v + " > " + v +".sorted")
  	      	timer=randint(1,6)
        	time.sleep(timer)
                os.system("mv "+v +".sorted " +v)
