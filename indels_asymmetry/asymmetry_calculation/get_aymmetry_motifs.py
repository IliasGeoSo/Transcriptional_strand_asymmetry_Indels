import re,os,sys,glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pdb
from scipy.stats import sem

def reader_file(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def reader_columns(path):
	Data=reader_file(path)
	Data2=[]
	DataR=[]
	for i in Data:
		DataR+=[int(i[-2])]
		Data2+=[int(i[-1])]
        return DataR,Data2

def reader_columns_boot(Data):
        Data2=[]
        DataR=[]
        for i in Data:
        	DataR+=[int(i[-2])]
                Data2+=[int(i[-1])]
        return DataR,Data2

def bootstrap_resampling(DataL):
        Bootstrap_DataL=[]
        size=len(DataL)
        import random    
        for i in range(size):
                chosen=random.choice(DataL)
                Bootstrap_DataL+=[chosen]
        return Bootstrap_DataL
		
# cancers to use
cancers = [ 'Bone','Bladder','CNS','Ovary', 'Breast', 'Skin', 'Kidney', 'Prostate', 'Esophagus', 'Pancreas', 'Lymphoid', 'Lung', 'Colorectal', 'Liver', 'Stomach', 'Biliary', 'Uterus','Head_neck']

colors = ["r","b","g","y","c"]

# output asymmetry scores
datafile=open("data_polyT_cancers","w")

# list of motifs for which to calculate mutational strand asymmetry
motifs = [["TT","TTT","TTTT","TTTTT","TTTTTT","TTTTTTT","TTTTTTTT","TTTTTTTTT","TTTTTTTTTT"]]
Ratios=[]
Sim_Ratios=[]
for cancer in cancers:
	for dk,motif_family in enumerate(motifs):
		plus_plus_minus_minusL=[]
		plus_plus_minus_minus_ML=[]
		plus_minus_minus_plusL=[]
		plus_minus_minus_plus_ML=[]

		All_line=[]
		scores1=[]
		scores2=[]
		for dkd, motif in enumerate(motif_family):
			semsL=[]
			# plus plus - non-template 
			files1_m=glob.glob("genes.detail.+.strand_plus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_plus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_plus_strand")[0]
			# minus minus - non-template
			files2_m=glob.glob("genes.detail.-.strand_minus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_minus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_minus_strand")[0]
			# minus plus - template
			files3_m=glob.glob("genes.detail.-.strand_plus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_minus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_minus_strand")[0]
			# plus minus - template
			files4_m=glob.glob("genes.detail.+.strand_minus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_plus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_plus_strand")[0]
			plus_plus_minus_minus = reader_columns(files1_m)[0]+reader_columns(files2_m)[0]
			plus_plus_minus_minus_M = reader_columns(files1_m)[1]+reader_columns(files2_m)[1]
			plus_minus_minus_plus = reader_columns(files4_m)[0]+reader_columns(files3_m)[0]
			plus_minus_minus_plus_M = reader_columns(files4_m)[1]+reader_columns(files3_m)[1]
	
			plus_plus_minus_minusL+=[sum(plus_plus_minus_minus)*len(motif)]
			plus_plus_minus_minus_ML+=[sum(plus_plus_minus_minus_M)]
			plus_minus_minus_plusL+=[sum(plus_minus_minus_plus)*len(motif)]
			plus_minus_minus_plus_ML+=[sum(plus_minus_minus_plus_M)]

		# Calculate asymmetry score
		Ratio_Real_part1 = sum(plus_plus_minus_minus_ML)/float(sum(plus_plus_minus_minusL))
		Ratio_Real_part2 = sum(plus_minus_minus_plus_ML) / float(sum(plus_minus_minus_plusL))
		Ratio_Real = Ratio_Real_part1 / (Ratio_Real_part1+Ratio_Real_part2)
		Ratios+=[Ratio_Real]

		# write asymmetry to output file
		datafile.write("polyT 2-10nt"+'\t'+str(sum(plus_plus_minus_minusL))+'\t'+str(sum(plus_minus_minus_plusL))+'\t'+str(sum(plus_plus_minus_minus_ML))+'\t'+str(sum(plus_minus_minus_plus_ML))+'\t'+str(Ratio_Real_part1)+'\t'+str(Ratio_Real_part2)+'\t'+str(Ratio_Real)+'\t'+cancer+'\n')


        #Bootstrapping here
	All_Ratios=[]
	for mk in range(1000):
		Ratios_Boot=[]
		for dk,motif_family in enumerate(motifs):
 	                plus_plus_minus_minusLB=[]
        	        plus_plus_minus_minus_MLB=[]
        	        plus_minus_minus_plusLB=[]
        	        plus_minus_minus_plus_MLB=[]
			for dkd, motif in enumerate(motif_family):
				# plus plus - non-template
        	                files1_m=glob.glob("genes.detail.+.strand_plus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_plus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_plus_strand")[0]
				# minus minus - non-template
        	                files2_m=glob.glob("genes.detail.-.strand_minus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_minus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_minus_strand")[0]
				# minus plus - template
        	                files3_m=glob.glob("genes.detail.-.strand_plus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_minus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_minus_strand")[0]
				# plus minus - template
        	                files4_m=glob.glob("genes.detail.+.strand_minus_genome.intersect_c."+motif+"_Coordinates_Genome.bed.genes_plus_strand.intersect_u."+cancer+"*indels*.intersect_u."+motif+"_Coordinates_Genome.bed.genes_plus_strand")[0]

        	                plus_plus_minus_minus = reader_columns(files1_m)[0]+reader_columns(files2_m)[0]
        	                plus_plus_minus_minus_M = reader_columns(files1_m)[1]+reader_columns(files2_m)[1]
        	                plus_minus_minus_plus = reader_columns(files4_m)[0]+reader_columns(files3_m)[0]
        	                plus_minus_minus_plus_M = reader_columns(files4_m)[1]+reader_columns(files3_m)[1]
	
                		plus_plus_minus_minus_Boot,plus_plus_minus_minus_M_Boot=reader_columns_boot(bootstrap_resampling(reader_file(files1_m)+reader_file(files2_m)))
                	        plus_minus_minus_plus_Boot,plus_minus_minus_plus_M_Boot=reader_columns_boot(bootstrap_resampling(reader_file(files3_m)+reader_file(files4_m)))
	
 		                plus_plus_minus_minusLB+=[sum(plus_plus_minus_minus_Boot)*len(motif)]
                	        plus_plus_minus_minus_MLB+=[sum(plus_plus_minus_minus_M_Boot)]
                       		plus_minus_minus_plusLB+=[sum(plus_minus_minus_plus_Boot)*len(motif)]
                        	plus_minus_minus_plus_MLB+=[sum(plus_minus_minus_plus_M_Boot)]
			
			# Calculate asymmetry score here	
		        Ratio_Real_part1_b = sum(plus_plus_minus_minus_MLB)/float(sum(plus_plus_minus_minusLB))
                	Ratio_Real_part2_b = sum(plus_minus_minus_plus_MLB) / float(sum(plus_minus_minus_plusLB))
                	Ratio_Real_b = Ratio_Real_part1_b / (Ratio_Real_part1_b+Ratio_Real_part2_b)
			Ratios_Boot+=[Ratio_Real_b]
		All_Ratios+=[Ratios_Boot]

	All_Ratios_S=np.array(All_Ratios).T.tolist()
	Sim_Ratios+=[All_Ratios_S]
datafile.close()


# Generate plot of asymmetry scores
matplotlib.rcParams['figure.facecolor'] = 'white'
ax = plt.subplot(111)
Medians = [np.median(k) for k in Ratios]
# sort 
Ratios_N, cancers = zip(*sorted(zip(Ratios, cancers)))
Ratios, Sim_Ratios = zip(*sorted(zip(Ratios, Sim_Ratios)))
plt.errorbar(Ratios, range(1,len(Ratios)+1,1), xerr=[np.std(mna) for mna in Sim_Ratios], linestyle='None',fmt='--o',color="green")#label="indels")
plt.axvline(x=0.5, hold=None,color="black")

# change cancer name from Head_neck to Head / Neck for display
cancers2=[]
for i in cancers:
	if "Head" in i:
		cancers2+=["Head / Neck"]	
	else:	
		cancers2+=[i]
cancers=cancers2

plt.yticks(range(1,len(cancers)+1,1),cancers)
plt.ylim(0,len(cancers)+1)
plt.title("polyT tracts")
plt.xlim(0.40,0.60)
plt.xlabel("T <-----------------------------------------  strand bias  -----------------------------------------> NT")
plt.grid()
ax.get_yaxis().tick_left()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("motifs_T2_10_all_cancers2.png")
plt.close()
