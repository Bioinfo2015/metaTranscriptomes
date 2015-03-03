#!/usr/bin/env python
##Import Modules 
from sys import exit
from sys import argv

orthos = argv[1]
dnds = argv[2]
outFile = argv[3]

orthoList = []
with open(orthos, 'r') as infile:
    for line in infile:
        orthoList.append(line.strip('\n'))
dnList = []
dsList = []
index = -1
orthoList2 = []
with open(dnds, 'r') as infile:
    for line in infile:
        index += 1
        line = line.strip('\n').split()
        try:
            dn = line[10]
            ds = line[13]
        except IndexError:
            dn = line[10]
            ds = line[12].split('=')[1]
        dnList.append(dn)
        dsList.append(ds)
        orthoList2.append(orthoList[index])

if len(orthoList) != len(dnList):
    exit("\nError, unequal Number of Orthologs and dN/dS values\n")

with open(outFile, 'w') as out:
    out.write("CNIDB\tdN\tdS")
    for i in range(len(orthoList)):
        out.write("\n{}\t{}\t{}".format(orthoList[i], dnList[i], dsList[i]))
        
        
    



                
            
            

