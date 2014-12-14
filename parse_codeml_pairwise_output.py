#!/usr/bin/env python
##parse_codeml_pairwise_output.py
##written 6/26/14 by Groves Dixon
ProgramName = 'parse_codeml_pairwise_output.py'
LastUpdated = '6/26/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Parses a list of codeml output files that were generated using pair-wise
dN/dS estimation (runmode -2). Pairs are set up against one base species
(set as spp1) and all other species (a list file)

'''

AdditionalProgramInfo = '''
Additional Program Information:


'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-f', required = True, dest = 'files', nargs="+", help = 'A glob to the codeml output files (probably *.codeml)')
parser.add_argument('-spp1', required = True, dest = 'spp1', help = 'The search tag for species 1')
parser.add_argument('-sppList', required = True, dest = 'sppList', help = 'The List of species to pair with species 1')
parser.add_argument('-o', required = True, dest = 'out', help = 'The desired output file name')
args = parser.parse_args()

#Assign Arguments
FileList = args.files
Spp1 = args.spp1
SppListName = args.sppList
OutfileName = args.out
SppList = []
with open(SppListName, 'r') as infile:
    for line in infile:
        SppList.append(line.strip("\n"))


def read_files(FileList, Spp1, SppList):
    '''Function to reads through each file and parses out
    dN and dS estimates for the specified species pair.
    '''
    print "\nLooking for data in {} codeml output files...".format(len(FileList))
    geneList = []
    dNList = []
    dSList = []
    speciesList = []
    for species in SppList:
        if species == Spp1:
            continue
        for file in FileList:
            with open(file, 'r') as infile:
                hit = 0
                hitCount = 0 #this should never exceed 1
                for line in infile:
                    if hitCount > 1:
                        exit("Found more than one instance of pairing in a file. Something is wrong.")
                    if hit == 0:
                        ##look for your species pair
                        if "("+Spp1+")" in line:
                            if "("+species+")" in line:
                                if "..." in line:
                                    hit = 1
                                    continue
                    elif hit == 1:
                        if "dN/dS=" in line:
                            line = line.split()
                            try:
                                dn = line[10]
                                ds = line[13]
                            except IndexError: #occurs sometimes when dS is very large
                                print line
                                #the dn value is also sometimes so high it must be split differently
                                #this probably means its a bad alignment/ortholog call, but pasrse it anyway
                                try:
                                    dn = line[10]
                                    ds = line[12].split('=')[1] #split the large ds value assuming that dS is >= 10.0 but dN is not
                                except IndexError:
                                    dn = line[9].split('=')[1] #this means that the dN value was also >= 10.0, so grab it differently
                                    ds = line[11].split('=')[1] #dS is also in a different place because of the big dN, so grab it
                            geneName = file.strip(".codeml")
                            geneList.append(geneName)
                            dNList.append(dn)
                            dSList.append(ds)
                            speciesList.append(species)
                            hit = 0
                            hitCount += 1
                            # print geneName
                            # print species
                            # print dn
    return geneList, dNList, dSList, speciesList
                        
def output(OutfileName, geneList, dNList, dSList, speciesList):
    """Outputs the data into a table"""
    with open(OutfileName, 'w') as out:
        out.write("EST\tspecies\tdN\tdS")
        for i in range(len(geneList)):
            out.write("\n{}\t{}\t{}\t{}".format(geneList[i], speciesList[i], dNList[i], dSList[i]))         
                            

geneList, dNList, dSList, speciesList = read_files(FileList, Spp1, SppList)
output(OutfileName, geneList, dNList, dSList, speciesList)


#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


