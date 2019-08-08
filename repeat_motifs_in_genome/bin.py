import re,os,sys,glob
import numpy as np

#generate the binned parts of genes, 14 bins, 10 bins between TSS and TES per gene. 2 bins upstream and 2 bins downstream of the TSS 10kB each.

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		Data+=[i.strip().split('\t')]
	return Data

def size(path):
	DataL=reader(path)
	Sizes=[]
	for v in DataL:
		Sizes+=[int(v[2])-int(v[1])]
	return Sizes

# Genic bins (We use 10 bins here)
Data_plus=reader("genes.detail.+.strand.txt")
Data_minus=reader("genes.detail.-.strand.txt")
bins=10
for binned in range(1,bins+1):
	datafile1=open("gene_coordinates_bin_"+str(binned)+"_plus.txt","w")
	datafile2=open("gene_coordinates_bin_"+str(binned)+"_minus.txt","w")
	for v in Data_plus:
		chrom = v[0]
		start = int(v[1])
		end = int(v[2])
		size = int(v[2])-int(v[1])
		size_bin = size /bins
		start_bin = start + (binned-1)*size_bin
		end_bin = start + (binned-1)*size_bin + size_bin
		datafile1.write(chrom+'\t'+str(start_bin)+'\t'+str(end_bin)+'\n')
	
	for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
		end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
		datafile2.write(chrom+'\t'+str(end_bin)+'\t'+str(start_bin)+'\n')
datafile1.close()
datafile2.close()

# Upstream 10kB
datafile1=open("gene_coordinates_upstream_10000nt_plus.txt","w")
datafile2=open("gene_coordinates_upstream_10000nt_minus.txt","w")
for v in Data_plus:
                chrom = v[0]
                start = int(v[1])
                end = int(v[2])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start + (binned-1)*size_bin
                end_bin = start + (binned-1)*size_bin + size_bin
		if start-10000>0:
        	        datafile1.write(chrom+'\t'+str(start-10000)+'\t'+str(start)+'\n')

for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
                end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
                datafile2.write(chrom+'\t'+str(end)+'\t'+str(end+10000)+'\n')
datafile1.close()
datafile2.close()

# Upstream 10kB to 20kB
datafile1=open("gene_coordinates_upstream_10000nt_20000nt_plus.txt","w")
datafile2=open("gene_coordinates_upstream_10000nt_20000nt_minus.txt","w")
for v in Data_plus:
                chrom = v[0]
                start = int(v[1])
                end = int(v[2])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start + (binned-1)*size_bin
                end_bin = start + (binned-1)*size_bin + size_bin
		if start-20000>0:
        	        datafile1.write(chrom+'\t'+str(start-20000)+'\t'+str(start-10000)+'\n')

for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
                end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
                datafile2.write(chrom+'\t'+str(end+10000)+'\t'+str(end+20000)+'\n')
datafile1.close()
datafile2.close()

# Downstream 10kB
datafile1=open("gene_coordinates_downstream_10000nt_plus.txt","w")
datafile2=open("gene_coordinates_downstream_10000nt_minus.txt","w")
for v in Data_plus:
                chrom = v[0]
                start = int(v[1])
                end = int(v[2])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start + (binned-1)*size_bin
                end_bin = start + (binned-1)*size_bin + size_bin
                datafile1.write(chrom+'\t'+str(end)+'\t'+str(end+10000)+'\n')

for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
                end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
		if start-10000>0:
        	        datafile2.write(chrom+'\t'+str(start-10000)+'\t'+str(start)+'\n')
datafile1.close()
datafile2.close()

# Downstream 10kB to 20kB
datafile1=open("gene_coordinates_downstream_10000nt_20000nt_plus.txt","w")
datafile2=open("gene_coordinates_downstream_10000nt_20000nt_minus.txt","w")
for v in Data_plus:
                chrom = v[0]
                start = int(v[1])
                end = int(v[2])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start + (binned-1)*size_bin
                end_bin = start + (binned-1)*size_bin + size_bin
                datafile1.write(chrom+'\t'+str(end+10000)+'\t'+str(end+20000)+'\n')

for v in Data_minus:
                chrom = v[0]
                start = int(v[2])
                end = int(v[1])
                size = int(v[2])-int(v[1])
                size_bin = size /bins
                start_bin = start - (binned-1)*size_bin
                end_bin = start - (binned-1)*size_bin - size_bin
		if start-20000>0:
        	        datafile2.write(chrom+'\t'+str(start-20000)+'\t'+str(start-10000)+'\n')
datafile1.close()
datafile2.close()

