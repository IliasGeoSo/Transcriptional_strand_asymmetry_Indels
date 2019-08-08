import re,os,sys,glob

datafile=open("polyNs.fa","w")
size=range(2,11)
for s in size:
	datafile.write(">"+str("A"*s)+'\n')
	datafile.write(s*"A"+'\n')

        datafile.write(">"+str("T"*s)+'\n')
        datafile.write(s*"T"+'\n')

        datafile.write(">"+str("C"*s)+'\n')
        datafile.write(s*"C"+'\n')
	
        datafile.write(">"+str("G"*s)+'\n')
        datafile.write(s*"G"+'\n')
datafile.close()
