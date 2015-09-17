#!/usr/bin/env python
#parse_codonW_effectiveNumberCodons.py
##Import Modules 
from sys import exit
from sys import argv

file = argv[1]
lineNumber = 0
print "contig\tNc"
with open(file, 'r') as infile:
    for line in infile:
        lineNumber += 1
        if lineNumber == 1:
            continue
        line = line.strip("\n").split()
        line[0] = line[0].split("_")[0]
        try:
            float(line[1])
        except ValueError:
            continue
        print("\t".join(line))
    
    

        
        
                
    
    

