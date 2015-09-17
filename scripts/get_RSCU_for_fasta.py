#!/usr/bin/env python
##get_RSCU_for_fasta.py
##written 9/8/15 by Groves Dixon
ProgramName = 'get_RSCU_for_fasta.py'
LastUpdated = '9/8/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This program reads through a fasta file of coding sequences.
It divides the sequence into codons assuming the reading frame begins with the 
first nucleotide.
It outputs a table with each seuqnce ID linked to the RSCU (relative synonymous codon usage)
for each codon and the species name (66 column table).
'''

AdditionalProgramInfo = '''
Additional Program Information:
Note the fasta must be RNA.

This script has been tested against output from dnaSP.

Ambiguous codons that include any uncertain nucleotide
such as N, W or Y are ignored in the analysis.
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
parser.add_argument('-i',    required = True, dest = 'input', help = 'The the input file name')
parser.add_argument('-o',    required = True, dest = 'output', help = 'The the desired output file name')
# parser.add_argument('-spp',  required = True, dest = 'spp',   help = 'The name of the species the data are for')
args = parser.parse_args()


#Usage
#get_RSCU_for_fasta -i sequences.fasta -o codon_RSCUs.txt


#Assign Arguments
InfileName  = args.input
OutfileName = args.output
debug = False
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
codonList = geneticCode.keys()   #global variable of all codons
codonList.sort()

#use the genetic code dictionary above to build 
#a new one to store counts in
def create_count_dict():
    """Function to make a clean dictionary 
    to keep counts of codons in"""
    codonList = geneticCode.keys()
    codonList.sort()
    countsDict = {}
    for i in geneticCode.keys():
        countsDict[i] = 0
    return countsDict


#we also need the reverse dictionary, where each list of synonymous 
#codons is paired with its amino acid
def build_inverse_code(geneticCode):
    """Function to revserse the genetic code entered above
    so that amino acids key to lists of synonymous codons
    that code for them"""
    aaList = []
    revCode = {}
    for i in codonList:
        aa = geneticCode[i]
        if aa not in aaList:
            aaList.append(aa)
            revCode[aa] = [i]
        else:
            revCode[aa].append(i)
    return aaList, revCode


def get_counts(InfileName):
    '''Function to read the fasta file and 
    count the codon use
    '''
    geneList = []
    dataDict = {}
    totalAmbiguous = 0 #keep count of ambiguous codons that have 'N's in them
    totalCodons = 0
    fasSeqs = SeqIO.parse(open(InfileName), 'fasta')
    for seqRecord in fasSeqs:
        seq = str(seqRecord.seq)
        seqID = seqRecord.id
        geneList.append(seqID)
        codonIndices = range(0, len(seq), 3)
        seqCodons = []                           #list of the codons in this sequence
        dataDict[seqID] = create_count_dict()    #set up empty counts dictionary for this gene
        for i in codonIndices:
            totalCodons += 1
            codon = seq[i:(i+3)].upper()
            #skip ambiguous codons that have an N in them
            if "N" in codon:
                totalAmbiguous += 1
                continue
            #skip ambiguous codons that have other characters
            try:
                geneticCode[codon]
            except KeyError:
                print "Warning, codon {} is not in standard genetic code and will be ignored".format(codon)
                totalAmbiguous += 1
                continue
            seqCodons.append(codon)
            dataDict[seqID][codon] += 1
    print "\n{} out of {} Codons were ambiguous and not scored in the RSCU calculations".format(totalAmbiguous, totalCodons)
    return geneList, dataDict


def get_aa_totals(revCode, aaList, geneticCode, codonList, codonCountDict):
    """Function to get the total number of times an amino 
    acid appears in a sequence by totaling the counts of its
    synonymous codons"""
    aaTotals = {}
    for aa in aaList:
        aaTotals[aa] = 0
    for codon in codonList:
        codonCounts = codonCountDict[codon]
        aa = geneticCode[codon]
        aaTotals[aa] += codonCounts
    return(aaTotals)
        
        

def calulate_rscu(aaList, revCode, geneList, dataDict):
    """Go through each gene and convert the counts data to RSCU"""
    rscuDict = {}
    for gene in geneList:
        rscuDict[gene] = {}
        codonCountDict = dataDict[gene]
        aaTotals = get_aa_totals(revCode, aaList, geneticCode, codonList, codonCountDict)
        for codon in codonList:
            codonCount = codonCountDict[codon]
            aa = geneticCode[codon]
            numberSynonymous = float(len(revCode[aa]))
            aaCount = float(aaTotals[aa])
            if aaCount == 0:
                rscuDict[gene][codon] = "NA"
                continue
            expected = aaCount / numberSynonymous
            #rather than cound missing codons as zero, all codons get a minimum count of 1
            # if codonCount == 0:
            #     codonCount += 1
            rscu = float(codonCount) / expected
            rscuDict[gene][codon] = str(np.round(rscu, decimals = 2))
    return rscuDict
            
def organize_codons(geneticCode, codonList, revCode, aaList):
    aaList.sort()
    sortedCodonList = []
    sortedAAs = []
    for aa in aaList:
        sCodons = revCode[aa]
        for sc in sCodons:
            sortedCodonList.append(sc)
            sortedAAs.append(aa)
    return sortedCodonList, sortedAAs

def output(OutfileName, rscuDict, geneticCode, revCode, geneList, sortedCodonList, sortedAAs):
    """Function to output the data as a table"""
    headerList = ['contig']
    header2List = ['contig']
    for i in sortedCodonList:
        headerList.append(i)
    for i in sortedAAs:
        header2List.append(i)
    if debug == True:
        print
        print
        print len(headerList)
        print headerList
        print len(header2List)
        print header2List
    with open(OutfileName, 'w') as out:
        header = "\t".join(headerList)
        header2 = "\t".join(header2List)
        out.write(header)
        # out.write("\n" + header2)
        for gene in geneList:
            rscuValues = []
            for codon in sortedCodonList:
                rscuValues.append(rscuDict[gene][codon])
            dataList = [gene] + rscuValues
            out.write("\n" + "\t".join(dataList))
            

aaList, revCode = build_inverse_code(geneticCode)
geneList, dataDict = get_counts(InfileName)
rscuDict = calulate_rscu(aaList, revCode, geneList, dataDict)
sortedCodonList, sortedAAs = organize_codons(geneticCode, codonList, revCode, aaList)
output(OutfileName, rscuDict, geneticCode, revCode, geneList, sortedCodonList, sortedAAs)



#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


