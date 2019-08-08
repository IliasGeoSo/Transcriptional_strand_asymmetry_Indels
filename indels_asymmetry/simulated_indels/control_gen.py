import re,os,sys
import glob
from random import shuffle
import random

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def in_here(coordinates):
        DataL=reader2("hg19.chrom.sizes")
        counts=0
        for i in DataL:
                if i[0]==coordinates[0]:
                        if int(coordinates[1])>0 and int(coordinates[1])<int(i[1]):
                                counts+=1
        return counts

def read_file(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Seqs=[]
        for i in data[1:]:
                Seqs+=[i.strip()]
        Seq="".join(Seqs)
        return Seq.upper()

def times_found_in_list2(Data,DataL):
        """"
        DataL has start end in pos 1 and pos 2
	Queries if two coordinates overlap
        """
        times_found=0
        for i in DataL:
		#testing if same chromosome
                if Data[0]==i[0]:
			#testing if one coordinate overlaps the other coordinate 
                        if (int(Data[1])>=int(i[1]) and int(Data[1])<=int(i[2]))  or (int(Data[2])>=int(i[1]) and int(Data[2])<=int(i[2])) or (int(i[1])>=int(Data[1]) and int(i[1])<=int(Data[2])) or  (int(i[2])>=int(Data[1]) and int(i[2])<=int(Data[2])):
                                times_found+=1
        return times_found

#reader non-mappable regions for mutation callling. Provide path
BlackList=reader("/nfs/compgen-04/team218/short_term_igs/masked_sites/All_non_mappable")
# Saves each path to each cancer type. All files are bed formatted. Provide path to indels.
files= glob.glob("/nfs/compgen-04/team218/ilias/indels_pancan_motif_analysis_ins_dels_by_tissue/samples/*.indels")

for one in files:	
	name = one.split("/")[-1].split(".indels")[0]
	if not os.path.exists(name):
		os.makedirs(name)
	Data=reader(one)
	# Number of jobs for parallel computing
        parts=4500
        total=len(Data)
        part=(int(sys.argv[1])-1)
        step=int(total*float(1/float(parts)))
	datafile=open(name+"/"+name+".simulations_indels_"+str(part),"w")
        DataLToSimulateL=Data[part*step:(part+1)*step]
        if part==parts-1:
                DataLToSimulateL=Data[part*step:]
        if part>10001 or part<0:
                break
	DataL=DataLToSimulateL

	for seq in DataL:
		#region to generate the simulated indel in
		window=2500
		chromosome_to_search =seq[0]
		position_center = (int(seq[1])+int(seq[2]))/2
	
		new_position = random.choice(range(position_center-window,position_center+window))
		indel_size = abs(int(seq[2])-int(seq[1]))
		PositionsControl=[chromosome_to_search,new_position,new_position+indel_size]

		while times_found_in_list2(PositionsControl,BlackList)>0 or PositionsControl==[] or len(PositionsControl)==0 or in_here(PositionsControl)==0:
                        chromosome_to_search = seq[0]
			new_position = random.choice(range(position_center-window,position_center+window))
			PositionsControl=[chromosome_to_search,new_position,new_position+indel_size]
		else:
			datafile.write(str(PositionsControl[0])+'\t'+str(int(PositionsControl[1]))+'\t'+str(int(PositionsControl[2]))+'\n')
	datafile.close()
				
