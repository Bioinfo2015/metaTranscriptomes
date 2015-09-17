#!/usr/bin/env python
##get_CAI_from_fasta.py
##written 9/15/15 by Groves Dixon
ProgramName = 'get_CAI_from_fasta.py'
LastUpdated = '9/15/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script caluclates the codon adaptation index (CAI) for a set of sequences
based on an input set of relative adaptiveness values (wi) for each codon.
Relative adaptiveness values can be generated using calculate_relative_adaptiveness.py.

Outputs a table of sequence IDs from the fasta file as rows linked with their
repective CAI values.
'''

AdditionalProgramInfo = '''
Additional Program Information:

See example of relative adaptiveness table at bottom of script.
'''

##Import Modules 
import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from scipy import stats
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input fasta file you want CAI values for')
parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file.')
parser.add_argument('-w', required = True, dest = 'wi', help = 'A table of relative adaptiveness values (wi) for each codon')
parser.add_argument('-c', required = False, default = 3, dest = 'col', help = 'The column number where the wi values are in the relative adaptiveness table. Codons are assumed to be in the first column. Default = 3rd column.')
args = parser.parse_args()

#Assign Arguments
debug = False
InfileName = args.input
OutfileName = args.out
wiFile = args.wi
wiColumn = int(args.col) - 1
nucleotideList = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']


def read_wi_table(wiFile):
    """Reads in the relative adaptiveness file"""
    wiDict = {}
    lineNumber = 0
    with open(wiFile, 'r') as infile:
        for line in infile:
            lineNumber += 1
            if lineNumber == 1: continue #skip header
            line = line.strip("\n").split("\t")
            codon = line[0]
            wi = line[wiColumn]
            wiDict[codon] = wi
    return wiDict
    

def get_cai(InfileName, wiDict):
    '''Function that reads through a set of coding sequences
    in fasta format and gets the geometric mean of the wi
    values for each codon in the sequence.
    '''
    numSeqs = 0
    resultsDict = {}
    seqIdList = []
    caiList = []
    fasSeqs = SeqIO.parse(open(InfileName), 'fasta')
    for seq in fasSeqs:
        numSeqs += 1
        seqString = str(seq.seq)
        codonCount = len(seqString)/3
        indexSeries = range(0, codonCount * 3, 3)
        codonList = []
        wiList = []
        for i in indexSeries:
            #some codons will have ambiguous nucleotides, so skip those
            good = 1
            codon = seqString[i:i+3].upper()
            for nucleotide in codon:
                if nucleotide in nucleotideList:
                    continue
                else:
                    good = 0
            if good == 1:
                codonList.append(codon)
        for codon in codonList:
            wiList.append(float(wiDict[codon]))
        cai = stats.gmean(wiList)
        #record the results in two lists for output
        caiList.append(cai)
        seqIdList.append(seq.id)
        if debug == True:
            print "\n------------"
            print "gene = {}".format(seq.id)
            print "codonList:"
            print codonList
            print "wiList:"
            print wiList
            print "CAI = {}".format(cai)
    return seqIdList, caiList
                
def output(seqIdList, caiList):
    """Output the results as a table"""
    with open(OutfileName, 'w') as out:
        out.write('contig\tcai')
        for i in range(len(seqIdList)):
            outString = "{}\t{}".format(seqIdList[i], caiList[i])
            out.write('\n' + outString)
        
                



wiDict = read_wi_table(wiFile)
seqIdList, caiList = get_cai(InfileName, wiDict)
output(seqIdList, caiList)


#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


#example of input wi table
#these were caluclated from a set of the top 5% 
#highest expressed genes in an A.millepora RNAseq dataset
exampleRelativeAdaptivenessTable = """
Example Relative Adaptiveness Table:
codon	aa	wi	rscu
GCA	Ala	0.956834532374	1.33
GCC	Ala	0.575539568345	0.8
GCG	Ala	0.338129496403	0.47
GCU	Ala	1.0	1.39
AGA	Arg	1.0	1.98
AGG	Arg	0.570707070707	1.13
CGA	Arg	0.484848484848	0.96
CGC	Arg	0.338383838384	0.67
CGG	Arg	0.222222222222	0.44
CGU	Arg	0.414141414141	0.82
AAC	Asn	0.834862385321	0.91
AAU	Asn	1.0	1.09
GAC	Asp	0.652892561983	0.79
GAU	Asp	1.0	1.21
UGC	Cys	0.754385964912	0.86
UGU	Cys	1.0	1.14
CAA	Gln	1.0	1.06
CAG	Gln	0.88679245283	0.94
GAA	Glu	1.0	1.21
GAG	Glu	0.652892561983	0.79
GGA	Gly	1.0	1.56
GGC	Gly	0.519230769231	0.81
GGG	Gly	0.339743589744	0.53
GGU	Gly	0.711538461538	1.11
CAC	His	0.754385964912	0.86
CAU	His	1.0	1.14
AUA	Ile	0.478571428571	0.67
AUC	Ile	0.671428571429	0.94
AUU	Ile	1.0	1.4
CUA	Leu	0.344827586207	0.5
CUC	Leu	0.475862068966	0.69
CUG	Leu	0.793103448276	1.15
CUU	Leu	0.965517241379	1.4
UUA	Leu	0.565517241379	0.82
UUG	Leu	1.0	1.45
AAA	Lys	1.0	1.09
AAG	Lys	0.834862385321	0.91
AUG	Met	1.0	1.0
UUC	Phe	0.639344262295	0.78
UUU	Phe	1.0	1.22
CCA	Pro	1.0	1.63
CCC	Pro	0.39263803681	0.64
CCG	Pro	0.282208588957	0.46
CCU	Pro	0.779141104294	1.27
AGC	Ser	0.698529411765	0.95
AGU	Ser	0.882352941176	1.2
UCA	Ser	1.0	1.36
UCC	Ser	0.558823529412	0.76
UCG	Ser	0.389705882353	0.53
UCU	Ser	0.882352941176	1.2
ACA	Thr	1.0	1.53
ACC	Thr	0.490196078431	0.75
ACG	Thr	0.346405228758	0.53
ACU	Thr	0.777777777778	1.19
UGG	Trp	1.0	1.0
UAC	Tyr	0.960784313725	0.98
UAU	Tyr	1.0	1.02
GUA	Val	0.467153284672	0.64
GUC	Val	0.569343065693	0.78
GUG	Val	0.890510948905	1.22
GUU	Val	1.0	1.37
"""

