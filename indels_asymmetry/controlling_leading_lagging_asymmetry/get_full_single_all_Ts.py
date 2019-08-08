import re,os,sys,glob,pdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np
import pandas as pd
import pdb
from scipy.stats import sem

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		Data+=[i.strip().split('\t')]
	return len(Data)

def reader2(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def reader3(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
	Data2=[]
	DataR=[]
	for i in Data:
		DataR+=[int(i[-2])]
		Data2+=[int(i[-1])]
        return DataR,Data2

def reader3_b(Data):
        Data2=[]
        DataR=[]
        for i in Data:
		DataR+=[int(i[-2])]
                Data2+=[int(i[-1])]
        return DataR,Data2


def bootstrap_resampling(DataL):
        size=len(DataL)
        import random    
        for i in range(size):
                chosen=random.choice(DataL)
                Bootstrap_DataL+=[chosen]
        return Bootstrap_DataL
		
ScoresL=[]
cancers = [ 'CNS','Ovary', 'Breast', 'Skin', 'Kidney', 'Prostate', 'Esophagus', 'Pancreas', 'Lymphoid', 'Lung', 'Colorectal', 'Liver', 'Stomach', 'Biliary', 'Uterus','Head_neck']
colors = ["r","b","g","y","c"]
motifs = [["TT","TTT","TTTT","TTTTT","TTTTTT","TTTTTTT","TTTTTTTT","TTTTTTTTT","TTTTTTTTTT"]]
stats =["MCF7_RepliStrand.leading","MCF7_RepliStrand.lagging"]
RatiosL=[]
Sim_RatiosL=[]
for stat in stats:
	Ratios=[]
	Sim_Ratios=[]
	for cancer in cancers:
		plus_plus_minus_minusL=[]
                plus_plus_minus_minus_ML=[]
                plus_minus_minus_plusL=[]
                plus_minus_minus_plus_ML=[]
		print cancer
		for dk,motif_family in enumerate(motifs):

			All_line=[]
			scores1=[]
			scores2=[]
			for dkd, motif in enumerate(motif_family):
				semsL=[]
				#mutations
				files1_m=glob.glob("genes.detail.+.strand_plus_genome.i_c.pgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.pgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat)[0]
				files2_m=glob.glob("genes.detail.-.strand_minus_genome.i_c.mgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.mgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat)[0]
				files3_m=glob.glob("genes.detail.+.strand_minus_genome.i_c.mgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.mgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat)[0]
				files4_m=glob.glob("genes.detail.-.strand_plus_genome.i_c.pgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.pgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat)[0]
				plus_plus_minus_minus = reader3(files1_m)[0]+reader3(files2_m)[0]
				plus_plus_minus_minus_M = reader3(files1_m)[1]+reader3(files2_m)[1]
				plus_minus_minus_plus = reader3(files4_m)[0]+reader3(files3_m)[0]
				plus_minus_minus_plus_M = reader3(files4_m)[1]+reader3(files3_m)[1]
	
				plus_plus_minus_minusL+=[sum(plus_plus_minus_minus)]
				plus_plus_minus_minus_ML+=[sum(plus_plus_minus_minus_M)]
				plus_minus_minus_plusL+=[sum(plus_minus_minus_plus)]
				plus_minus_minus_plus_ML+=[sum(plus_minus_minus_plus_M)]
	

		Ratio_Real_part1 = sum(plus_plus_minus_minus_ML)/float(sum(plus_plus_minus_minusL))
		Ratio_Real_part2 = sum(plus_minus_minus_plus_ML) / float(sum(plus_minus_minus_plusL))
		Ratio_Real = Ratio_Real_part1 / float((Ratio_Real_part1+Ratio_Real_part2))
		Ratios+=[Ratio_Real]
		print sum(plus_plus_minus_minus_ML),sum(plus_minus_minus_plus_ML)
		print Ratios
        	#Bootstrapping
		All_Ratios=[]
		for mk in range(200):
			Ratios_Boot=[]
			for dk,motif_family in enumerate(motifs):
 	        	        plus_plus_minus_minusLB=[]
        		        plus_plus_minus_minus_MLB=[]
        		        plus_minus_minus_plusLB=[]
        		        plus_minus_minus_plus_MLB=[]
				for dkd, motif in enumerate(motif_family):
                                	files1_m=glob.glob("genes.detail.+.strand_plus_genome.i_c.pgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.pgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat)[0]
                                	files2_m=glob.glob("genes.detail.-.strand_minus_genome.i_c.mgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.mgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat)[0]
                                	files3_m=glob.glob("genes.detail.+.strand_minus_genome.i_c.mgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.mgeno."+motif+"_Coordinates_Genome.bed.genes_plus_strand.i_u."+stat)[0]
                                	files4_m=glob.glob("genes.detail.-.strand_plus_genome.i_c.pgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat+".i_u."+cancer+".indels.i_u.pgeno."+motif+"_Coordinates_Genome.bed.genes_minus_strand.i_u."+stat)[0]

        		                plus_plus_minus_minus = reader3(files1_m)[0]+reader3(files2_m)[0]
        	        	        plus_plus_minus_minus_M = reader3(files1_m)[1]+reader3(files2_m)[1]
        	        	        plus_minus_minus_plus = reader3(files4_m)[0]+reader3(files3_m)[0]
        	        	        plus_minus_minus_plus_M = reader3(files4_m)[1]+reader3(files3_m)[1]
		
        	        		plus_plus_minus_minus_Boot,plus_plus_minus_minus_M_Boot=reader3_b(bootstrap_resampling(reader2(files1_m)+reader2(files2_m)))
        	        	        plus_minus_minus_plus_Boot,plus_minus_minus_plus_M_Boot=reader3_b(bootstrap_resampling(reader2(files3_m)+reader2(files4_m)))
		
 			                plus_plus_minus_minusLB+=[sum(plus_plus_minus_minus_Boot)]
                		        plus_plus_minus_minus_MLB+=[sum(plus_plus_minus_minus_M_Boot)]
                	       		plus_minus_minus_plusLB+=[sum(plus_minus_minus_plus_Boot)]
                        		plus_minus_minus_plus_MLB+=[sum(plus_minus_minus_plus_M_Boot)]
				

		        	Ratio_Real_part1_b = sum(plus_plus_minus_minus_MLB)/float(sum(plus_plus_minus_minusLB))
                		Ratio_Real_part2_b = sum(plus_minus_minus_plus_MLB) / float(sum(plus_minus_minus_plusLB))
                		Ratio_Real_b = Ratio_Real_part1_b / (Ratio_Real_part1_b+Ratio_Real_part2_b)
				Ratios_Boot+=[Ratio_Real_b]
			All_Ratios+=[Ratios_Boot]

		All_Ratios_S=np.array(All_Ratios).T.tolist()
		Sim_Ratios+=[All_Ratios_S]
		print '\n'

	RatiosL+=[Ratios]
	Sim_RatiosL+=[Sim_Ratios]

plt.rcParams['axes.facecolor'] = 'white'
ax = plt.subplot(111)
#Leading strand
Medians = [np.median(k) for k in RatiosL[0]]
Ratios_N, cancers = zip(*sorted(zip(RatiosL[0], cancers)))
RatiosL_S, Sim_RatiosL_S = zip(*sorted(zip(RatiosL[0], Sim_RatiosL[0])))
plt.errorbar(RatiosL_S, range(1,len(RatiosL[0])+1,1), xerr=[np.std(mna) for mna in Sim_RatiosL_S], linestyle='None',fmt='--o',color="lightblue",label="leading")
#Lagging strand
Medians = [np.median(k) for k in RatiosL[1]]
Ratios, Sim_Ratios = zip(*sorted(zip(RatiosL[1], Sim_RatiosL[1])))
plt.errorbar(Ratios, range(1,len(RatiosL[0])+1,1), xerr=[np.std(mna) for mna in Sim_Ratios], linestyle='None',fmt='--o',color="lightpink",label="lagging")
plt.axvline(x=0.5,hold=None,color="black")
plt.yticks(range(1,len(cancers)+1,1),cancers,fontsize=14)
plt.xlabel("T <-------------------------------  strand bias  -------------------------------> NT",fontsize=14)
plt.ylim(0,len(cancers)+0.5)
plt.xlim(0.4,0.6)
plt.grid()
ax.get_yaxis().tick_left()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.title("polyTs",fontsize=16)
plt.tight_layout()
plt.legend(loc="lower right",frameon=False,numpoints=1)
plt.savefig("motifs_T2_10_all_cancers_leading_lagging.png")
plt.close()
