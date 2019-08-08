import re,os,sys,glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random
from scipy.stats import sem
import pandas as pd
import pdb

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		Data+=[i.strip().split('\t')]
	return Data

def nucs(path):
	Data=reader(path)
	Nucs=0
	for v in Data:
		Nucs+=int(v[2])-int(v[1])
	return Nucs

def boot(FirstL,SecondL):
        from random import choice
        # take each consecutive element from both lists and see if that ratio is expected
        Ratios_errs=[]
        for k in range(len(FirstL)):
                Ratios=[]
                for v in range(100):
                        Occs=[]
                        element_first = int(FirstL[k])
                        element_second = int(SecondL[k])
    
                        Picks = [ 0 for dj in range(element_first)]+[ 1 for dk in range(element_second)]
    
                        for d in range(element_first+element_second):
                                Occs+=[choice(Picks)]
                        if Occs.count(1)!=0:
                                Ratios+=[Occs.count(0)/float(Occs.count(1))]
                        else:
                                pass
                Ratios_errs+=[sem(Ratios)]
        return Ratios_errs

motifs=["TT","TTT","TTTT","TTTTT","TTTTTT","TTTTTTT","TTTTTTTT","TTTTTTTTT","TTTTTTTTTT","AA","AAA","AAAA","AAAAA","AAAAAA","AAAAAAA","AAAAAAAA","AAAAAAAAA","AAAAAAAAAA","GG","GGG","GGGG","GGGGG","GGGGGG","GGGGGGG","GGGGGGGG","CC","CCC","CCCC","CCCCC","CCCCCC","CCCCCCC","CCCCCCCC"]
parts = ["gene_coordinates_upstream_10000nt_20000nt_minus","gene_coordinates_upstream_10000nt_minus","gene_coordinates_bin_1_minus","gene_coordinates_bin_2_minus","gene_coordinates_bin_3_minus","gene_coordinates_bin_4_minus","gene_coordinates_bin_5_minus","gene_coordinates_bin_6_minus","gene_coordinates_bin_7_minus","gene_coordinates_bin_8_minus","gene_coordinates_bin_9_minus","gene_coordinates_bin_10_minus","gene_coordinates_downstream_10000nt_minus","gene_coordinates_downstream_10000nt_20000nt_minus"]

ExpectedR=[]
score_non_templateL=[]
score_templateL=[]
for mot in motifs:
	PerMotif=[]
	PerMotifR=[]
	PerMotif_tL=[]
	PerMotif_ntL=[]
	for v in parts:
		# read plus plus orientation genic file -> non-template
		filed_plus_plus = reader(glob.glob(mot+"_Coordinates_Genome.bed_genome_plus.intersect_u.*"+v.split("minus")[0]+"plus.txt")[0])
		# read minus minus orientation genic file -> non-template
		filed_minus_minus = reader(glob.glob(mot+"_Coordinates_Genome.bed_genome_minus.intersect_u.*"+v.split("minus")[0]+"minus.txt")[0])
		# read minus plus orientation genic file -> template
		filed_minus_plus = reader(glob.glob(mot+"_Coordinates_Genome.bed_genome_minus.intersect_u.*"+v.split("minus")[0]+"plus.txt")[0])
		# read plus minus orientation genic file -> template
		filed_plus_minus = reader(glob.glob(mot+"_Coordinates_Genome.bed_genome_plus.intersect_u.*"+v.split("minus")[0]+"minus.txt")[0])
		score_non_template=len(filed_plus_plus)+len(filed_minus_minus)
		score_template=len(filed_minus_plus)+len(filed_plus_minus)
		PerMotif+=[float(score_non_template)/float(score_template)]
		PerMotifR+=[float(score_non_template)/float(score_template)]
		PerMotif_tL+=[float(score_template)]
		PerMotif_ntL+=[float(score_non_template)]
        ExpectedR+=[PerMotifR]
        score_non_templateL+=[PerMotif_ntL]
        score_templateL+=[PerMotif_tL]

# heatmap enrichment of 
df = pd.DataFrame(ExpectedR,index=motifs,columns=["upstream 20kB","upstream 10kB"]+["bin " +str(dv) for dv in range(1,11,1)]+["downstream 10kB","downstream 20kB"])
g=sns.heatmap(df,vmin=0,vmax=2,cmap='coolwarm')
for item in g.get_xticklabels():
    item.set_rotation(60)
plt.xlabel("Gene length",fontsize=12)
plt.tight_layout()
plt.savefig("heatmap_background_occs_motifs_bias.png",transparent=True)
plt.close()

# alter naming
motif_names=[]
for i in motifs:
        motif_names+=[i[0]+"${"+str(len(i))+"}$"]

df = pd.DataFrame(ExpectedR,index=motif_names,columns=["upstream 20kB","upstream 10kB"]+["bin " +str(dv) for dv in range(1,11,1)]+["downstream 10kB","downstream 20kB"])
g=sns.heatmap(df,vmin=0,vmax=2,cmap='coolwarm')
for item in g.get_xticklabels():
    item.set_rotation(90)
plt.ylabel("polyN length",fontsize=12)
plt.xlabel("Gene length",fontsize=12)
plt.tight_layout()
plt.savefig("heatmap_background_occs_motifs_bias2.png",transparent=True)
plt.close()

min2=nucs("../polyN/gene_coordinates_upstream_10000nt_20000nt.txt")
min1=nucs("../polyN/gene_coordinates_upstream_10000nt.txt")
get1=nucs("../polyN/gene_coordinates_bin_1.txt")
get2=nucs("../polyN/gene_coordinates_bin_2.txt")
get3=nucs("../polyN/gene_coordinates_bin_3.txt")
get4=nucs("../polyN/gene_coordinates_bin_4.txt")
get5=nucs("../polyN/gene_coordinates_bin_5.txt")
get6=nucs("../polyN/gene_coordinates_bin_6.txt")                   
get7=nucs("../polyN/gene_coordinates_bin_7.txt")                   
get8=nucs("../polyN/gene_coordinates_bin_8.txt")
get9=nucs("../polyN/gene_coordinates_bin_9.txt")
get10=nucs("../polyN/gene_coordinates_bin_10.txt")
max1=nucs("../polyN/gene_coordinates_downstream_10000nt.txt")
max2=nucs("../polyN/gene_coordinates_downstream_10000nt_20000nt.txt")
backsL=[min2,min1,get1,get2,get3,get4,get5,get6,get7,get8,get9,get10,max1,max2]
All_BDens=[]
for motif in motifs:
	Total_Occs=0
	for indexs,part in enumerate(parts):
		filed_plus_plus = glob.glob(motif+"_Coordinates_Genome.bed_genome_plus.intersect_u.*"+part.split("minus")[0]+"*plus.txt")[0]
                filed_minus_minus = glob.glob(motif+"_Coordinates_Genome.bed_genome_minus.intersect_u.*"+part.split("minus")[0]+"*minus.txt")[0]
                filed_minus_plus = glob.glob(motif+"_Coordinates_Genome.bed_genome_minus.intersect_u.*"+part.split("minus")[0]+"*plus.txt")[0]
                filed_plus_minus = glob.glob(motif+"_Coordinates_Genome.bed_genome_plus.intersect_u.*"+part.split("minus")[0]+"*minus.txt")[0]
		occsN=0
		occsN+=len(reader(filed_plus_plus))+len(reader(filed_minus_minus))+len(reader(filed_minus_plus))+len(reader(filed_plus_minus))
		Total_Occs+=occsN
	All_BDens+=[Total_Occs/float(sum(backsL))]

All_Dens=[]
for dnd,motif in enumerate(motifs):
        Dens=[]
        for indexs,part in enumerate(parts):
                filed_plus_plus = glob.glob(motif+"_Coordinates_Genome.bed_genome_plus.intersect_u.*"+part.split("minus")[0]+"*plus.txt")[0]
                filed_minus_minus = glob.glob(motif+"_Coordinates_Genome.bed_genome_minus.intersect_u.*"+part.split("minus")[0]+"*minus.txt")[0]
                filed_minus_plus = glob.glob(motif+"_Coordinates_Genome.bed_genome_minus.intersect_u.*"+part.split("minus")[0]+"*plus.txt")[0]
                filed_plus_minus = glob.glob(motif+"_Coordinates_Genome.bed_genome_plus.intersect_u.*"+part.split("minus")[0]+"*minus.txt")[0]
                occsN=0
		occsN+=len(reader(filed_plus_plus))+len(reader(filed_minus_minus))+len(reader(filed_minus_plus))+len(reader(filed_plus_minus))
                if "upstream_10000nt_20000nt" in part:
			back = backsL[0]
                elif "upstream_10000nt" in part and "20000" not in part:
                        back = backsL[1]
		elif "bin_1_" in part:	
			back = backsL[2]
		elif "bin_2_" in part:
			back = backsL[3]
		elif "bin_3_" in part:
			back = backsL[4]
		elif "bin_4_" in part:
			back = backsL[5]
		elif "bin_5_" in part:
			back = backsL[6]
		elif "bin_6_" in part:
			back = backsL[7]
		elif "bin_7_" in part:
			back = backsL[8]
		elif "bin_8_" in part:
			back = backsL[9]
		elif "bin_9_" in part:
			back = backsL[10]
		elif "bin_10_" in part:	
			back = backsL[11]
		elif "downstream_10000nt" in part and "20000" not in part:
			back = backsL[12]
		elif "downstream_10000nt_20000nt" in part:
			back = backsL[13]
		Dens+=[(occsN/float(back))/float(All_BDens[dnd])]
	All_Dens+=[Dens]

motifs_J=["AA","AAA","AAAA","AAAAA","AAAAAA","AAAAAAA","AAAAAAAA","AAAAAAAAA","AAAAAAAAAA"]
motifs_J2=["GG","GGG","GGGG","GGGGG","GGGGGG","GGGGGGG","GGGGGGGG","GGGGGGGG"]
scoresL=[]
for mot in motifs_J:
        filed_plus_plus = glob.glob(mot+"_Coordinates_Genome.bed_genome_plus.intersect_u.genes.detail.+.strand.txt")[0]
        filed_minus_minus = glob.glob(mot+"_Coordinates_Genome.bed_genome_minus.intersect_u.genes.detail.-.strand.txt")[0]
        filed_minus_plus = glob.glob(mot+"_Coordinates_Genome.bed_genome_minus.intersect_u.genes.detail.+.strand.txt")[0]
        filed_plus_minus = glob.glob(mot+"_Coordinates_Genome.bed_genome_plus.intersect_u.genes.detail.-.strand.txt")[0]

	score_template=len(reader(filed_plus_minus))+len(reader(filed_minus_plus))
	score_non_template=len(reader(filed_plus_plus))+len(reader(filed_minus_minus))	

	scoresL+=[float(score_non_template)/float(score_template)]
	

motif_names=[]
for i in motifs:
        motif_names+=[i[0]+"${"+str(len(i))+"}$"]

df = pd.DataFrame(All_Dens,index=motif_names,columns=["upstream 20kB","upstream 10kB"]+["bin " +str(dv) for dv in range(1,11,1)]+["downstream 10kB","downstream 20kB"])
datafile=open("all_out","w")
for k in All_Dens:
	datafile.write('\t'.join([str(x) for x in k])+'\n')
datafile.close()

df.to_csv("output_heatmap_input_df", sep='\t')
g=sns.heatmap(df,vmin=0.5,vmax=1.5,cmap="RdBu")
for item in g.get_xticklabels():
    item.set_rotation(90,fontsize=14)
for item in g.get_yticklabels():
	item.set_rotation(fontsize=14)
plt.xlabel("Gene body",fontsize=14)
plt.ylabel("polyN length",fontsize=14)
plt.tight_layout()
plt.savefig("heatmap_background_occs_motifs_Enrichments.png")
plt.close()

df.to_csv("output_heatmap_input_df", sep='\t')



