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
hilo refers to the high and low ends of the principal axis
produced by the codonW correspondance analysis. The genes
at the high end are assumed to be highly expressed and lowly
expressed. The origentation of the principal axis is arbitrary,
to the high end could be the negative or positve end of the axis.
Which end is assigned as the high end is based on the effective 
number of codons, which should be lower for the genes showing greater
codon bias. Hense, high refers to the 5 percent of genes at the pole of 
principal axis that has the lower Nc. The low fraction refers to 
the opposite 5%.
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
parser.add_argument('-hilo', required = True, dest = 'hilo',  help = 'The desired output name for the hilo data')
parser.add_argument('-cai',  required = True, dest = 'cai',   help = 'The desired name for the cai data')
parser.add_argument('-spp',  required = True, dest = 'spp',   help = 'The name of the species the data are for')
args = parser.parse_args()

#Assign Arguments
debug = 0         #set to one to get more printed information as script runs
InfileName  = args.input
hiloOutName = args.hilo
caiOutName  = args.cai
species = args.spp

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

def split_usage(uLine):
    """The usage expressions in the summary file are not tab separated.
    This function deals with these phases individually and returns a list
    of the four pieces of data. Used in function parse_hilo_line"""
    uLine = uLine.split()
    hiUse = uLine[0]
    hiCount = uLine[1][1:-1]
    lowUse = uLine[2]
    lowCount = uLine[3][1:-1]
    return [hiUse, hiCount, lowUse, lowCount]
    

def prep_codon(codon):
    """Performs some opperations on the codon given in summary table
    of hilo data"""
    codon = codon.replace(" ", "")
    codon = codon.strip("@")
    optimal = '0'
    if "*" in codon:
        optimal = '1'
        codon = codon.strip("*")
    return codon, optimal


def cpg_codon(codon):
    cg = '0'
    gc = '0'
    if "CG" in codon:
        cg = '1'
    if "GC" in codon:
        gc = '1'
    return [cg, gc]
    
def parse_hilo_line(line):
    codon1, optimal1 = prep_codon(line[1])
    codon2, optimal2 = prep_codon(line[4])
    cg1, gc1 = cpg_codon(codon1)
    cg2, gc2 = cpg_codon(codon2)
    aa1 = geneticCode[codon1]
    aa2 = geneticCode[codon2]
    usage1 = line[2]
    usage2 = line[5]
    u1 = split_usage(usage1)
    u2 = split_usage(usage2)
    result1 = [codon1, aa1] + u1 + [optimal1] + [cg1, gc1]
    result2 = [codon2, aa2] + u2 + [optimal2] + [cg2, gc2]
    if debug == 1:
        print line      #original line
        print result1   #extacted data for first codon
        print result2   #extracted data for second codon
    return result1, result2
        
def get_hilo(InfileName):
    '''Function to read in a file as a list of lists
    '''
    print "\nGathering Hilo Data from Summary file..."
    lineNumber = 0
    record = 0
    resultsList = []
    codonCountsDict = {}
    with open(InfileName, 'r') as infile:
        for line in infile: 
            if record == 0:
                if "extremes of axis 1" in line:
                    record = 1
                    continue
            else:
                if "Number of codons in high bias dataset" in line:
                    hiCounts = line.split()[-1]
                    codonCountsDict['hiCounts'] = hiCounts
                    continue
                if "Number of codons in low  bias dataset" in line:
                    lowCounts = line.split()[-1]
                    codonCountsDict['lowCounts'] = lowCounts
                    print "\nFinished gathering hilo data"
                    break
                else:
                    line = line.strip("\n").split("\t")
                    if len(line) <= 1:
                        continue
                    else:
                        result1, result2 = parse_hilo_line(line)
                        if result1[1] != "STOP":
                            resultsList.append(result1)
                        if result2[1] != "STOP":
                            resultsList.append(result2)
    return resultsList, codonCountsDict

def output_hilo_data(hiloOutName, hiloData, codonCountsDict):
    """Output the hilo data extracted by get_hilo()"""
    with open(hiloOutName, 'w') as out:
        #header is not included so the files can be concatenated
        # out.write("Codon\tAA\tHiUse\tHiCount\tlowUse\tlowCount\toptimal\tCpG.codon\tGpCcodon\ttotalHiCodons\ttotalLowCodons\tspecies")
        for line in hiloData:
            data = line + [codonCountsDict['hiCounts'], codonCountsDict['lowCounts']] + [species]
            print data
            out.write("\n{}".format("\t".join(data)))


def get_cai(InfileName):
    """Function to gather the cai data for codons"""
    print "\nGathering CAI Data from Summary file..."
    lineNumber = 0
    record = 0
    resultsList = []
    codonCountsDict = {}
    with open(InfileName, 'r') as infile:
        for line in infile:
            if record == 0:
                if "Cod AA    Xi	Wi" in line:
                    record = 1
                    continue
                else: continue
            else:
                line = line.strip("\n").split("\t")
                set1 = line[0].split()
                set2 = line[1].split()
                codon1 = set1[0]
                codon2 = set2[0]
                cg1 = cpg_codon(codon1)
                cg2 = cpg_codon(codon2)
                set1 += cg1
                set2 += cg2
                resultsList.append(set1)
                resultsList.append(set2)
    print "\nFinished Gathering cai data"
    return resultsList
                
def output_cai(caiOutName, resulsList):
    with open(caiOutName, 'w') as out:
        #header "codon\tAA\tXi\tWi\tCG_codon\tGCcodon\tspecies"
        #heder is not included to allow for easy concatenating
        for line in resulsList:
            out.write("\n{}\t{}".format("\t".join(line), species))


hiloData, codonCountsDict = get_hilo(InfileName)
output_hilo_data(hiloOutName, hiloData, codonCountsDict)
caiData = get_cai(InfileName)
output_cai(caiOutName, caiData)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


