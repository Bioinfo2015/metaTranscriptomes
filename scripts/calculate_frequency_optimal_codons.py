#!/usr/bin/env python
##calculate_frequency_optimal_codons.py
##written 2/23/15 by Groves Dixon
ProgramName = 'calculate_frequency_optimal_codons.py'
LastUpdated = '2/23/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script caluclates the frequency of optimal codons (Fop) for each sequence in a fasta.
It requires that you have a table of optimal codons for each amino acid (see AdditionalProgramInfo)
Input fasta is assumed to be a coding sequence in frame so that the first nucleotide in each sequence
is the first position nucleotide in the first codon.
'''

AdditionalProgramInfo = '''
Additional Program Information:
Codons with ambiguous bases are skipped.
Example Optimal codon table:
AA	optimal
F	UUU
L	UUG
I	AUU
V	GUU
S	UCU
P	CCA
T	ACA
A	GCU
Y	UAU
H	CAU
Q	CAA
N	AAU
K	AAA
D	GAU
E	GAA
C	UGU
R	AGA
S	AGU
G	GGA
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
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input file')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-optimal', required = True, dest = 'optimal', help = 'A subset of sequence names to pull and concatenate')
args = parser.parse_args()

#Assign Arguments
InfileName = args.input
OutfileName = args.out
OptimalCodonFile = args.optimal
nucleotideList = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']

def read_optimal_table(OptimalCodonFile):
    """Reads in the optimal codon file and produces a list of optimal codons
    and a dictionary linking each single letter amino acid abbreviation to its
    optimal codon"""
    optiDict = {}
    optimalList = []
    with open(OptimalCodonFile, 'r') as infile:
        for line in infile:
            line = line.strip("\n").split()
            AA = line[0]
            optimal = line[1]
            optiDict[AA] = optimal
            optimalList.append(optimal)
    return optiDict, optimalList
    

def get_Fop(InfileName, optimalList):
    '''Function read through the fasta using biopython and calculate Fop for each sequence.
    Each sequence is split into its constituent codons. Then each codon is scored as either 
    an optimal codon or not based on the given table. The ratio of optimal codons to total 
    codons for each sequence gives its Fop.
    '''
    numSeqs = 0
    resultsDict = {}
    seqIdList = []
    fasSeqs = SeqIO.parse(open(InfileName), 'fasta')
    for seq in fasSeqs:
        optimalCount = 0
        suboptimalCount = 0
        ambiguousCodonCount = 0
        numSeqs += 1
        seqString = str(seq.seq)
        codonCount = len(seqString)/3
        indexSeries = range(0, codonCount * 3, 3)
        codonList = []
        for i in indexSeries:
            codonList.append(seqString[i:i+3])
        # print "\n\n"
        # print "length = {}".format(len(seqString))
        # print codonCount
        # print indexSeries
        # print seqString
        # print codonList
        for codon in codonList:
            good = 1
            for nucleotide in codon:
                if nucleotide in nucleotideList:
                    continue
                else:
                    good = 0
            if good == 1:
                if codon in optimalList:
                    optimalCount += 1
                else:
                    suboptimalCount += 1
            else:
                ambiguousCodonCount += 1
        Fop = float(optimalCount) / float(codonCount)
        results = [str(codonCount), str(optimalCount), str(suboptimalCount), str(ambiguousCodonCount), str(Fop)]
        resultsDict[seq.id] = results
        seqIdList.append(seq.id)
    return seqIdList, resultsDict
                
def output(OutfileName, seqIdList, resultsDict):
    """Output the results as a table"""
    with open(OutfileName, 'w') as out:
        out.write('EST\ttotal\toptimal\tsuboptimal\tambiguous\tFop')
        for i in seqIdList:
            results = i + '\t' + '\t'.join(resultsDict[i])
            out.write('\n' + results)
        
                



optiDict, optimalList = read_optimal_table(OptimalCodonFile)
seqIdList, resultsDict = get_Fop(InfileName, optimalList)
output(OutfileName, seqIdList, resultsDict)


#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


