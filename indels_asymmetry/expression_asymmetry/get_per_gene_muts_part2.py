import re,os,sys,glob

files1 = glob.glob("../../../plus/*Coordinates_Genome.bed")
files2 = glob.glob("../../../minus/*Coordinates_Genome.bed")
#provide real path
path_to_mutations="path"
files=["Breast_cell_of_origin_RNA_Seq_q1_minus_genes.txt","Breast_cell_of_origin_RNA_Seq_q2_minus_genes.txt","Breast_cell_of_origin_RNA_Seq_q3_minus_genes.txt","Breast_cell_of_origin_RNA_Seq_q1_plus_genes.txt","Breast_cell_of_origin_RNA_Seq_q2_plus_genes.txt","Breast_cell_of_origin_RNA_Seq_q3_plus_genes.txt"]

#Below these outputs were generated from the previous script
Gen=[]
for j in files1:
	for k in mutations:
		Gen+=[k.split("/")[-1]+".intersect_u."+j.split("/")[-1]+".plus_genome"]

for j in files2:
        for k in mutations:
                Gen+=[k.split("/")[-1]+".intersect_u."+j.split("/")[-1]+".minus_genome"]


#Generate a file with the genes and a column of the polyN motifs and the polyNs overlalping the mutations
# 1. for plus strand
for i in files:
	for j in files1:
		os.system("bedtools intersect -a " + i + " -b " + j + " -c -sorted > " + i.split("/")[-1]+".intersect_c."+j.split("/")[-1])
		for k in Gen:
			motif = k.split("intersect_u.")[1].split("_")[0]
			motif_used = j.split("/")[-1].split("_")[0]
			if motif==motif_used and "plus_genome" in k:
				os.system("bedtools intersect -a " + i.split("/")[-1]+".intersect_c."+j.split("/")[-1]+" -b " + k + " -c -sorted > " +i.split("/")[-1]+".intersect_c."+j.split("/")[-1]+".intersect_c."+k.split("/")[-1])

# 2. for minus strand
for i in files:
        for j in files2:
                os.system("bedtools intersect -a " + i + " -b " + j + " -c -sorted > " + i.split("/")[-1]+".intersect_c."+j.split("/")[-1])
                for k in Gen:
                        motif = k.split("intersect_u.")[1].split("_")[0]
                        motif_used = j.split("/")[-1].split("_")[0]
                        if motif==motif_used and "minus_genome" in k:
                                os.system("bedtools intersect -a " + i.split("/")[-1]+".intersect_c."+j.split("/")[-1]+" -b " + k + " -c -sorted > " +i.split("/")[-1]+".intersect_c."+j.split("/")[-1]+".intersect_c."+k.split("/")[-1])

