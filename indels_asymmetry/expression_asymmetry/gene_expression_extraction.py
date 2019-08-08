import re,os,sys,glob
import pandas as pd
from pandas import DataFrame
import numpy as np

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()	
	Data=[]
	for i in data:
		Data+=[i.strip().split('\t')]
	return Data

def match(GeneID,CancerType,Ensembl):
        CoordinatesLL=[]
        for i in range(len(GeneID)):
                CoordinatesL=[]
                for j in Ensembl:
                        if j[0]==GeneID[i]:
                                CoordinatesL=["chr"+j[1],j[2],j[3],j[6],j[4],CancerType[i]]
                if CoordinatesL!=[]:
                        CoordinatesLL+=[CoordinatesL]
        return CoordinatesLL

Data=reader("57epigenomes.RPKM.pc")
Headers=Data[0]
Rest=Data[1:]
df= DataFrame(Rest,columns=Headers)

Ovarian=np.array(df["E097"])
Liver=np.array(df["E118"])
Pancreatic=np.array(df["E098"])
Esophagus=np.array(df["E079"])
Pedriatic_Brain_Tumour=np.array(df["E082"])
Gastric_cancer=np.array(df["E094"])
Malignant_Lymphoma=np.array(df["E116"])
Breast_cancer=np.array(df["E028"])
Stomach=np.array(df["E094"])
Colon=np.array(df["E106"])
Lung=np.array(df["E096"])
K562=np.array(df["E123"])
GM12878=np.array(df["E116"])
HelaS3=np.array(df["E117"])
Foreskin_keratinocyte = np.array(df["E057"])
Foreskin_fibroblast = np.array(df["E053"])

Genes=np.array(df["gene_id"])
Ensembl=reader("Ensembl_v65.Gencode_v10.ENSG.gene_info")

Coordinates_data_F=match(Genes,Stomach,Ensembl)
datafile=open("Stomach_cell_of_origin_RNA_Seq.txt","w")
for i in Coordinates_data_F:
        datafile.write(str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])+'\t'+str(i[5])+'\n')
datafile.close()
