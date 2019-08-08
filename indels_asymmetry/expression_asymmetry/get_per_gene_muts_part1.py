import re,os,sys,glob

#files=glob.glob("/nfs/compgen-04/team218/ilias/indels_pancan_motif_analysis_strand_bias_polyN2/general/All_chr_hg19_TSA/plus"
	
# provide path to polyN counts at plus and minus orientation here
files1 = glob.glob("../../../plus/*Coordinates_Genome.bed")
files2 = glob.glob("../../../minus/*Coordinates_Genome.bed")
#provide real path here
path_to_mutations="path"
mutations = glob.glob(path_to_mutations+"/Breast.indels_patients_DELETIONS")

files=["Breast_cell_of_origin_RNA_Seq_q1_minus_genes.txt","Breast_cell_of_origin_RNA_Seq_q2_minus_genes.txt","Breast_cell_of_origin_RNA_Seq_q3_minus_genes.txt","Breast_cell_of_origin_RNA_Seq_q1_plus_genes.txt","Breast_cell_of_origin_RNA_Seq_q2_plus_genes.txt","Breast_cell_of_origin_RNA_Seq_q3_plus_genes.txt"]
for j in files1:
	for k in mutations:
		os.system("bedtools intersect -a " + k  + " -b " + j + " -u -sorted > " + k.split("/")[-1]+".intersect_u."+j.split("/")[-1]+".plus_genome")

for j in files2:
        for k in mutations:
                os.system("bedtools intersect -a " + k  + " -b " + j + " -u -sorted > " + k.split("/")[-1]+".intersect_u."+j.split("/")[-1]+".minus_genome")


