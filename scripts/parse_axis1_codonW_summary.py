#!/usr/bin/env python
#parse_axis1_codonW_summary.py
##Import Modules 
from sys import exit
from sys import argv

summaryFile = argv[1]


dataList = []
lineNumber = 0
ready = 0
SET = 0
go = 0
oneBack = 0
twoBack = 0
current = 0
goSign = ['The position of each gene by axis', '(see also genes.coa)']
for line in summaryFile:
    lineNumber += 1
    line = line.strip('\n')
    if go == 1:
        if not line.strip():
            go = 0
            continue
        else:
            data = line.split()
            entry = data[0].split("_")
            gene = entry[1]
            number = entry[0]
            data[0] = gene
            result = [number, gene, data[1]]
            dataList.append(result)
    else:
        twoBack = oneBack
        oneBack = current
        current = line
        if [twoBack, oneBack] == goSign:
            print 'Found'
            print lineNumber
            exit()
        
        
                
    
    

