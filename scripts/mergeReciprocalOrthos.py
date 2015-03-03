#!/usr/bin/env python
##Import Modules 
from sys import exit
from sys import argv
import argparse

Description = """
Script to merge reciprocal orthologs from multiple species
and a single database. The sequence IDs from the database should be
in the first column of each file that is merged. The ortholog seq ID
from the particular species should be in the second column.
"""
AdditionalProgramInfo = '''
Additional Program Information:

To use this script, first get two fasta files
that you think contain orthologs.
Blast each fasta against the other using out format 5 
so that biopython can parse the blast output file.
Now you have the two fasta files, and two blast output
files to pass as arguments to this script.
'''


parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-f', required = True, dest = 'files', nargs = "+", help = 'A glob to the files you want to merge.')
parser.add_argument('-c', required = False, default = "none", dest = 'cut', help = 'Indicate the proportion of files that must have a reciprocal ortholog in order to output it. If a base ortholog does not have above that value it will be skipped')
parser.add_argument('-r', required = True, dest = 'ref', help = 'Indicate the species name for the reference. This should be different than any of the query species names. Example = AdigitiferaREF')
parser.add_argument('-o', required = False, default = "orthologs.txt", dest = 'out', help = 'The desired output name')
args = parser.parse_args()

FileList = args.files
Cutoff = args.cut
RefSpecies = args.ref
OutName = args.out
print "\nMerging Ortholog Tables for the following {} Files:".format(len(FileList))
for i in FileList:
    print i


dataDict = {} #place to store the putative ortholog pairs
fileCount = 0 #keep count of the files
for f in FileList:
    entryList = [] 
    fileCount += 1
    lineNumber = 0
    with open(f, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split('\t')
            if lineNumber == 1:
                header = line
                continue
            ortho1 = line[0]
            ortho2 = line[1]
            entryList.append(ortho1)
            try:
                dataDict[ortho1].append(ortho2)
            except KeyError:
                dataDict[ortho1] = [ortho2]

#filter the data dictionary based on the -c cutoff
filteredList = []
if Cutoff != "none":
    cutCount = 0
    passCount = 0
    print "\nFiltering Orthologs Based on a Cutoff of {}...".format(Cutoff)
    Cutoff = float(Cutoff)
    for i in entryList:
        orthoList = dataDict[i]
        # print orthoList
        trimList = []
        for x in orthoList:
            if x != "none":
                trimList.append(x)
        proportion = float(len(trimList))/float(len(orthoList))
        if proportion >= Cutoff:
            # print i
            passCount += 1
            filteredList.append(i)
        else:
            cutCount += 1
    print "{} Genes Passed and will be output to the results table".format(passCount)
    print "{} Genes Failed and will not be output".format(cutCount)
    entryList = filteredList
        
        


sppList = []
for i in FileList:
    sppList.append(i.split('_')[1]) ##set up header so that columns are just species names
Header = [RefSpecies] + sppList
print "\nOutputting merged results for the following Species:"
for i in sppList:
    print i
HeaderString = '\t'.join(Header)
with open(OutName, 'w') as out:
    out.write(HeaderString)
    for i in entryList:
        resultString = i +'\t' + '\t'.join(dataDict[i])
        out.write("\n" + resultString)



                
            
            

