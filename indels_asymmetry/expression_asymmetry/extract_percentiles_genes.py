import re,os,sys,glob
import numpy as np

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		Data+=[i.strip().split('\t')]
	return Data

cancer="Breast"
expressions=[]
DataL = reader(cancer+"_cell_of_origin_RNA_Seq.txt")
for i in DataL:
	expressions +=[float(i[-1])]

p1=np.percentile(expressions,33)
p2=np.percentile(expressions,66)

datafile11=open(cancer+"_cell_of_origin_RNA_Seq_q1_plus_genes.txt","w")
datafile21=open(cancer+"_cell_of_origin_RNA_Seq_q2_plus_genes.txt","w")
datafile31=open(cancer+"_cell_of_origin_RNA_Seq_q3_plus_genes.txt","w")

datafile12=open(cancer+"_cell_of_origin_RNA_Seq_q1_minus_genes.txt","w")
datafile22=open(cancer+"_cell_of_origin_RNA_Seq_q2_minus_genes.txt","w")
datafile32=open(cancer+"_cell_of_origin_RNA_Seq_q3_minus_genes.txt","w")



for i in DataL:
	# plus oriented genes
	if i[-2]=="1":
		if float(i[-1])<p1:
			datafile11.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\n')
		elif float(i[-1])<p2:
			datafile21.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\n')
		else:
			datafile31.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\n')
datafile11.close()
datafile21.close()
datafile31.close()

for i in DataL:
	# minus oriented genes
        if i[-2]=="-1":
                if float(i[-1])<p1:
                        datafile12.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\n')
                elif float(i[-1])<p2:
                        datafile22.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\n')
                else:
                        datafile32.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\n')
datafile12.close()
datafile22.close()
datafile32.close()
