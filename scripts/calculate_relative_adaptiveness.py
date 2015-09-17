#!/usr/bin/env python
##calculate_relative_adaptiveness.py
##written 9/15/15 by Groves Dixon
ProgramName = 'calculate_relative_adaptiveness.py'
LastUpdated = '9/15/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This program calculates the relative adaptiveness (w) of each codon
from the RSCU values from a concatentated set of highly expressed genes.

For each codon i that codes for an amino acid j, the relative
adaptiveness wij is equal to the ratio of the RSCU for that
codon to that of the most abundant synonymous codon:

wij = RSCUi / RSCUimax

(This is equal to the ratio of their counts)
'''

AdditionalProgramInfo = '''
Additional Program Information:
To use this first get a set of highly expressed genes.
Then extract and concatenate their coding sequences in to a single fasta.
Calculate RSCU values from the concatenated sequence using get_RSCU_for_fasta.py.
eg:
get_RSCU_for_fasta.py -i top5_hiExpressed_contigs_concatenated.fasta -o top5_RSCU.txt

Use the RSCU values as input for this script.
'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i',    required = True, dest = 'input', help = 'The the input file name with RSCUs for each codon from a set of highly expressed genes')
parser.add_argument('-o',    required = True, dest = 'output', help = 'The the desired output file name. Each codon will be grouped in a table with its amino acid and its relative adaptiveness value')
args = parser.parse_args()


#Assign Arguments
InfileName  = args.input
OutfileName = args.output
debug = True
# Species = args.spp

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

#set up a list of the codons  
def read_rscu(InfileName):
    """Funciton to read in the RSCU values.
    Output is a dictionary linking amino acids with 
    two parallel lists, their codons, and the codon's 
    RSCU values extracted from the input file."""
    lineNumber = 0
    with open(InfileName, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip("\n").split("\t")
            if lineNumber == 1:
                codonList = line[1:] #gather the codons from the input file from the first line
                aaList = []          #set up paralell list of amino acids
                for i in codonList:
                    aa = geneticCode[i]
                    aaList.append(aa)
                rscuDict = {}
                for i in aaList:
                    rscuDict[i] = [[],[]] #set up dictionary to store the set of codons and rscu values for each amino acid
            if lineNumber == 2:
                rscuList = line[1:]
                for i in range(len(rscuList)):
                    codon = codonList[i]
                    aa = aaList[i]
                    rscu = rscuList[i]
                    if rscu == 0:
                        rscu = 0.5  #this is to prevent w values of zero, which would reduce CAI to zero (see Behura and Severson 2013; Biological Reviews 88:49-61)
                    rscuDict[aa][0].append(codon)
                    rscuDict[aa][1].append(rscu)
    if debug == True:
        print "\nResults from gathering RSCU data:"
        aas = rscuDict.keys()
        aas.sort()
        for i in aas:
            print i
            print rscuDict[i]
    return rscuDict

def calculate_w(rscuDict, OutfileName):
    """Function to calculate relative
    adaptiveness from the assembled rscu data"""
    #organize the amino acids and alphebetize them
    with open(OutfileName, 'w') as out:
        header = "codon\taa\twi\trscu"
        out.write(header)
        aaList = rscuDict.keys()
        aaList.remove('STOP')
        aaList.sort()
        for i in aaList:
            rscus = rscuDict[i][1]
            codons = rscuDict[i][0]
            xmax = float(max(rscus))
            if len(codons) != len(rscus):
                exit("Error, don't have an rscu for each codon for {}".format(i))
            for c in range(len(codons)):
                codon = codons[c]
                ru = float(rscus[c])
                w = ru/xmax
                outString = "\n{}\t{}\t{}\t{}".format(codon, i, w, ru)
                
                
            
        
             
                
            
rscuDict = read_rscu(InfileName)
calculate_w(rscuDict, OutfileName)



#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


