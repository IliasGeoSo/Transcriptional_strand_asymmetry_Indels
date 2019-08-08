import re,os,sys,glob

# Replication leading / lagging compartments
MCF7s=["MCF7_RepliStrand.leading","MCF7_RepliStrand.lagging"]
# PolyN files for non-template and template analysis. The folder represents the genome direction used e.g. /plus/ the ending represents the gene direction e.g. genes_plus_strand
files1 = glob.glob("../../plus/*Coordinates_Genome.bed.genes_plus_strand")
files2 = glob.glob("../../minus/*Coordinates_Genome.bed.genes_minus_strand")
files3 = glob.glob("../../minus/*Coordinates_Genome.bed.genes_plus_strand")
files4 = glob.glob("../../plus/*Coordinates_Genome.bed.genes_minus_strand")
mutations = glob.glob("/nfs/compgen-04/team218/ilias/indels_pancan_motif_analysis_ins_dels_by_tissue/samples/*indels")
files1_=[]
for one in files1:
	name = one.split("/")[-1].split("_Coordinates_")[0]
	if len(name)>1:
		for k in MCF7s:
			os.system("bedtools intersect -a "+ one + " -b " + k + " -u > pgeno." + one.split("/")[-1]+".i_u."+k)
			files1_+=["pgeno."+one.split("/")[-1]+".i_u."+k]

files2_=[]
for one in files2:
	name = one.split("/")[-1].split("_Coordinates_")[0]
	if len(name)>1:
        	for k in MCF7s:
        	        os.system("bedtools intersect -a "+ one + " -b " + k + " -u > mgeno." + one.split("/")[-1]+".i_u."+k)
			files2_+=["mgeno."+one.split("/")[-1]+".i_u."+k]

files3_ =[]
for one in files3:
	name = one.split("/")[-1].split("_Coordinates_")[0]
	if len(name)>1:
        	for k in MCF7s:
                	os.system("bedtools intersect -a "+ one + " -b " + k + " -u > mgeno." + one.split("/")[-1]+".i_u."+k)
			files3_+=["mgeno."+one.split("/")[-1]+".i_u."+k]

files4_=[]
for one in files4:
	name = one.split("/")[-1].split("_Coordinates_")[0]
	if len(name)>1:
        	for k in MCF7s:
        	        os.system("bedtools intersect -a "+ one + " -b " + k + " -u > pgeno." + one.split("/")[-1]+".i_u."+k)
			files4_+=["pgeno."+one.split("/")[-1]+".i_u."+k]

files1=files1_
files2=files2_
files3=files3_
files4=files4_

files=["genes.detail.+.strand.txt","genes.detail.-.strand.txt"]
for j in files1:
	os.system("bedtools intersect "  + " -a genes.detail.+.strand.txt -b " + j + " -c -sorted > " +  "genes.detail.+.strand_plus_genome.i_c."+j.split("/")[-1])
	for k in mutations:
		os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > "+k.split("/")[-1]+".i_u."+j.split("/")[-1])
		os.system("bedtools intersect -a " + "genes.detail.+.strand_plus_genome.i_c."+j.split("/")[-1] + " -b  "+k.split("/")[-1]+".i_u."+j.split("/")[-1] + " -c > " + "genes.detail.+.strand_plus_genome.i_c."+j.split("/")[-1]+".i_u." + k.split("/")[-1]+".i_u."+j.split("/")[-1])

for j in files2:
	os.system("bedtools intersect " + " -a genes.detail.-.strand.txt -b " +j+ " -c -sorted > " +  "genes.detail.-.strand_minus_genome.i_c."+j.split("/")[-1])
        for k in mutations:
                os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > "+k.split("/")[-1]+".i_u."+j.split("/")[-1])
                os.system("bedtools intersect -a " + "genes.detail.-.strand_minus_genome.i_c."+j.split("/")[-1] + " -b  "+k.split("/")[-1]+".i_u."+j.split("/")[-1]+ " -c > " + "genes.detail.-.strand_minus_genome.i_c."+j.split("/")[-1]+".i_u." + k.split("/")[-1]+".i_u."+j.split("/")[-1])


for j in files3:
	os.system("bedtools intersect -a  genes.detail.+.strand.txt -b " + j + " -c -sorted > " +  "genes.detail.+.strand_minus_genome.i_c."+j.split("/")[-1])
        for k in mutations:
                os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > "+k.split("/")[-1]+".i_u."+j.split("/")[-1])
                os.system("bedtools intersect -a " + "genes.detail.+.strand_minus_genome.i_c."+j.split("/")[-1] + " -b  "+k.split("/")[-1]+".i_u."+j.split("/")[-1] + " -c > " + "genes.detail.+.strand_minus_genome.i_c."+j.split("/")[-1]+".i_u." + k.split("/")[-1]+".i_u."+j.split("/")[-1])


for j in files4:
	os.system("bedtools intersect -a genes.detail.-.strand.txt -b " + j + " -c -sorted > " +  "genes.detail.-.strand_plus_genome.i_c."+j.split("/")[-1])
        for k in mutations:
                os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > "+k.split("/")[-1]+".i_u."+j.split("/")[-1])
                os.system("bedtools intersect -a " + "genes.detail.-.strand_plus_genome.i_c."+j.split("/")[-1] + " -b  "+k.split("/")[-1]+".i_u."+j.split("/")[-1]+ " -c > " + "genes.detail.-.strand_plus_genome.i_c."+j.split("/")[-1]+".i_u." + k.split("/")[-1]+".i_u."+j.split("/")[-1])



