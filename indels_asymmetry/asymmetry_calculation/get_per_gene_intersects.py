import re,os,sys,glob

# path to files with coordinates of polyN motifs that were interesected with gene coordinates in repeat_motifs_in_genome/motif_finding/
files1 = glob.glob("../plus/*Coordinates_Genome.bed.genes_plus_strand")
files2 = glob.glob("../minus/*Coordinates_Genome.bed.genes_minus_strand")
files3 = glob.glob("../minus/*Coordinates_Genome.bed.genes_plus_strand")
files4 = glob.glob("../plus/*Coordinates_Genome.bed.genes_minus_strand")


#indel files separate for each tissue are provided as a path
path_to_indels ="path"
mutations = glob.glob(path_to_indels+"/*indels")
#sort indel files to perform bedtools commands
for i in mutations:
	os.system("sort -k1,1 -k2,2n " + i + " > " + i+".sorted")
	os.system("mv " + i + ".sorted " + i)

# perform bedtools functions from here on and save files of different orientations separate. E.g. plus_gene and plus_genome are non-template type
# plus plus orientation -> non-tempalte
for j in files1:
	os.system("bedtools intersect "  + " -a genes.detail.+.strand.txt -b " + j + " -c -sorted > per_gene/" +  "genes.detail.+.strand_plus_genome.intersect_c."+j.split("/")[-1])
	for k in mutations:
		os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1])
		os.system("bedtools intersect -a " + "per_gene/genes.detail.+.strand_plus_genome.intersect_c."+j.split("/")[-1] + " -b  per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1] + " -c > " + "per_gene/genes.detail.+.strand_plus_genome.intersect_c."+j.split("/")[-1]+".intersect_u." + k.split("/")[-1]+".intersect_u."+j.split("/")[-1])

# minus minus orientation -> non-tempalte
for j in files2:
	os.system("bedtools intersect " + " -a genes.detail.-.strand.txt -b " +j+ " -c -sorted > per_gene/" +  "genes.detail.-.strand_minus_genome.intersect_c."+j.split("/")[-1])
        for k in mutations:
                os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1])
                os.system("bedtools intersect -a " + "per_gene/genes.detail.-.strand_minus_genome.intersect_c."+j.split("/")[-1] + " -b  per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1]+ " -c > " + "per_gene/genes.detail.-.strand_minus_genome.intersect_c."+j.split("/")[-1]+".intersect_u." + k.split("/")[-1]+".intersect_u."+j.split("/")[-1])


# plus minus orientation -> tempalte
for j in files3:
	os.system("bedtools intersect -a  genes.detail.+.strand.txt -b " + j + " -c -sorted > per_gene/" +  "genes.detail.+.strand_minus_genome.intersect_c."+j.split("/")[-1])
        for k in mutations:
                os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1])
                os.system("bedtools intersect -a " + "per_gene/genes.detail.+.strand_minus_genome.intersect_c."+j.split("/")[-1] + " -b  per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1] + " -c > " + "per_gene/genes.detail.+.strand_minus_genome.intersect_c."+j.split("/")[-1]+".intersect_u." + k.split("/")[-1]+".intersect_u."+j.split("/")[-1])

# minus plus orientation - > template
for j in files4:
	os.system("bedtools intersect -a genes.detail.-.strand.txt -b " + j + " -c -sorted > per_gene/" +  "genes.detail.-.strand_plus_genome.intersect_c."+j.split("/")[-1])
        for k in mutations:
                os.system("bedtools intersect -a " + k + " -b " + j +" -u -sorted > per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1])
                os.system("bedtools intersect -a " + "per_gene/genes.detail.-.strand_plus_genome.intersect_c."+j.split("/")[-1] + " -b  per_gene/"+k.split("/")[-1]+".intersect_u."+j.split("/")[-1]+ " -c > " + "per_gene/genes.detail.-.strand_plus_genome.intersect_c."+j.split("/")[-1]+".intersect_u." + k.split("/")[-1]+".intersect_u."+j.split("/")[-1])



