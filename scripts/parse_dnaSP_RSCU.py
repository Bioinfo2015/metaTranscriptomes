#!/usr/bin/env python
##parse_codonW_hilo_and_cai.py
##written 9/1/15 by Groves Dixon
ProgramName = 'parse_codonW_hilo_and_cai.py'
LastUpdated = '9/1/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Extract the usage and codon adaptation indices for 
codons from a summary.coa file output from codonW.
'''

AdditionalProgramInfo = '''
Additional Program Information:
IMPORTANT. There is a bug in this script in that it cannot parse 
any sequence for which the number of uses for a particular codon is > 999.
So it will work for most sequences, but will break with large concatenated sequnces.
If that happens use 
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
parser.add_argument('-i',    required = True, dest = 'input', help = 'The the input file name')
parser.add_argument('-o',    required = True, dest = 'rscuOut', help = 'The the desired output file name')
# parser.add_argument('-optOut',    required = True, dest = 'optOut', help = 'The the desired optimal codon output file name')
parser.add_argument('-spp',  required = True, dest = 'spp',   help = 'The name of the species the data are for')
args = parser.parse_args()

#Assign Arguments
InfileName  = args.input
RscuOutName = args.rscuOut
# OptOutName = args.optOut
Species = args.spp

#set up genetic code
geneticCode = {"UUU":"Phe", "UUC":"Phe", "UUA":"Leu", "UUG":"Leu",
    "UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser",
    "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
    "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
    "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu",
    "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    "CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
    "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
    "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
    "AGU":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
    "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
    "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}

def get_rscu(InfileName):
    '''Function to read in a file as a list of lists
    '''
    print "\nGathering RSCU Data..."
    lineNumber = 0
    record = 0
    #set up lists to fill with the usage data
    codonList = []
    countList = []
    rscuList = []
    aaList = []
    codonIndices = [0, 3, 6, 9]  #locations of codon in each data line
    countIndices = [1, 4, 7, 10] #locations of the total counts for the codons
    rscuIndices = [2, 5, 8, 11]  #locations of the rscu values
    with open(InfileName, 'r') as infile:
        for line in infile: 
            if record == 0:
                if "RSCU and Codon Usage Table of" in line:
                    record = 1
                    continue
            else:
                line = line.strip("\n").split()
                if len(line) <= 1:
                    continue
                else:
                    print line
                    for i in codonIndices:
                        print i
                        codon = line[i].split("-")[0]
                        codonList.append(codon)
                        aaList.append(geneticCode[codon])
                    for i in rscuIndices:
                        ru = line[i][1:-1]
                        rscuList.append(ru)
                    for i in countIndices:
                        countList.append(line[i])
    return codonList, countList, rscuList, aaList



def output_rscu(RscuOutName, codonList, rscuList, aaList):
    """Output the hilo data extracted by get_hilo()"""
    with open(RscuOutName, 'w') as out:
        for i in range(len(codonList)):
            outString = "\n{}\t{}\t{}\t{}\t{}".format(codonList[i], countList[i], aaList[i], rscuList[i], Species)
            out.write(outString)


codonList, countList, rscuList, aaList = get_rscu(InfileName)
output_rscu(RscuOutName, codonList, rscuList, aaList)



#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


