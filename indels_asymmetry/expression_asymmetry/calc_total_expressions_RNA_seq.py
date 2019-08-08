import re,os,sys,glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		Data+=[i.strip().split('\t')]
	return Data

def dens(path):
	DataL=reader(path)
	total=0
	overlaps=0
	for i in DataL:
		total+=int(i[2])-int(i[1])
		overlaps+=int(i[-1])
	return overlaps,total

def dens2(path):
        DataL=reader(path)
        total=0
        overlaps=0
        for i in DataL:
                total+=int(i[2])-int(i[1])
                overlaps+=int(i[-2])
        return overlaps,total

motifs = ["GG","GGG","GGGG","GGGGG","GGGGGG","GGGGGGG","GGGGGGGG","GGGGGGGGG","GGGGGGGGGG"]
All_Sets=[]
cold=["r","g","b","y","k"]
cold = ["lightpink","hotpink","purple"]
lines = ["Lung","Liver","Breast","Esophagus","Ovarian","Stomach","Skin","Colon","K562","Pancreatic","CNS"]
cancers = ["Lung","Liver","Breast","Esophagus","Ovary","Stomach","Skin","Colorectal","Lymphoid","Pancreas","CNS"]
for mdm in range(len(lines)):
	cancer = cancers[mdm]
	line = lines[mdm]
	Set=[]
	for k in range(1,4):
		pp_mm_overlaps_All = 0
		pm_mp_overlaps_All = 0
		pp_mm_total_All = 0
		pm_mp_total_All = 0

		pp_mm_overlaps_b_All = 0 
	        pm_mp_overlaps_b_All = 0 
	        pp_mm_total_b_All = 0 
	        pm_mp_total_b_All = 0
		for motif in motifs:
			#path to files to use for strand asymmetry
			p_p=line+"_cell_of_origin_RNA_Seq_q"+str(k)+"_plus_genes.txt.intersect_c."+motif+"_Coordinates_Genome.bed.intersect_c."+cancer+".indels_patients_DELETIONS.intersect_u."+motif+"_Coordinates_Genome.bed.plus_genome"
			m_m=line+"_cell_of_origin_RNA_Seq_q"+str(k)+"_minus_genes.txt.intersect_c."+motif+"_Coordinates_Genome.bed.intersect_c."+cancer+".indels_patients_DELETIONS.intersect_u."+motif+"_Coordinates_Genome.bed.minus_genome"
			p_m=line+"_cell_of_origin_RNA_Seq_q"+str(k)+"_plus_genes.txt.intersect_c."+motif+"_Coordinates_Genome.bed.intersect_c."+cancer+".indels_patients_DELETIONS.intersect_u."+motif+"_Coordinates_Genome.bed.minus_genome"
			m_p=line+"_cell_of_origin_RNA_Seq_q"+str(k)+"_minus_genes.txt.intersect_c."+motif+"_Coordinates_Genome.bed.intersect_c."+cancer+".indels_patients_DELETIONS.intersect_u."+motif+"_Coordinates_Genome.bed.plus_genome"
	
			p_p_overlaps,p_p_total=dens(p_p)
			m_m_overlaps,m_m_total=dens(m_m)
			p_m_overlaps,p_m_total=dens(p_m)
			m_p_overlaps,m_p_total=dens(m_p)
	
			p_p_overlaps_b,p_p_total_b=dens2(p_p)
		       	m_m_overlaps_b,m_m_total_b=dens2(m_m)
        		p_m_overlaps_b,p_m_total_b=dens2(p_m)
	        	m_p_overlaps_b,m_p_total_b=dens2(m_p)
	
			pp_mm_overlaps = p_p_overlaps+m_m_overlaps
			pp_mm_total = p_p_total+m_m_total
			pm_mp_overlaps = p_m_overlaps+m_p_overlaps
			pm_mp_total = p_m_total+m_p_total
	
       			pp_mm_overlaps_b = p_p_overlaps_b+m_m_overlaps_b
       			pp_mm_total_b = p_p_total_b+m_m_total_b
        		pm_mp_overlaps_b = p_m_overlaps_b+m_p_overlaps_b
        		pm_mp_total_b = p_m_total_b+m_p_total_b

			pp_mm_overlaps_All+=pp_mm_overlaps
			pm_mp_overlaps_All+=pm_mp_overlaps
			pp_mm_total_All+=pp_mm_total
			pm_mp_total_All+=pm_mp_total
	
	        	pp_mm_overlaps_b_All+=pp_mm_overlaps_b
        		pm_mp_overlaps_b_All+=pm_mp_overlaps_b
      	  		pp_mm_total_b_All+=pp_mm_total_b
        		pm_mp_total_b_All+=pm_mp_total_b

		Dens_pp_mm = pp_mm_overlaps_All
		Dens_pm_mp = pm_mp_overlaps_All
        	Dens_pp_mm_b = pp_mm_overlaps_b_All
        	Dens_pm_mp_b = pm_mp_overlaps_b_All

		Ratio1 = Dens_pp_mm / float(Dens_pp_mm_b)
		Ratio2 = Dens_pm_mp / float(Dens_pm_mp_b)

		Set+=[Ratio1/float(Ratio1+Ratio2)]

	All_Sets+=[Set]

exps = ["Low","Medium","High"]
for mdm in range(len(All_Sets)):
        plt.scatter([mdm+1,mdm+1,mdm+1],All_Sets[mdm],color=cold[:3])
	plt.scatter([mdm+1,mdm+1,mdm+1],All_Sets[mdm],color=cold)

exps=exps[::-1]
cold=cold[::-1]
for k in range(3):
        plt.scatter([],[],color=cold[k],label=exps[k])
    
ax = plt.subplot(111)
plt.ylim(0.35,0.65)
plt.axhline(0.5)
plt.ylabel("Non-template strand asymmetry")
plt.grid()
plt.title("Deletions")
plt.xticks(range(1,len(lines)+1),cancers,rotation=90)
ax.get_yaxis().tick_left()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.legend(loc='lower right',scatterpoints=1,frameon=False,title="Expression levels")
plt.savefig("expression_dels_cancers_polyG.png")
plt.close()




	
	

	

