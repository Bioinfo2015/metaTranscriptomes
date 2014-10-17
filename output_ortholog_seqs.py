#!/usr/bin/env python
##ortholog_outputter.py
##written 9/16/14 by Groves Dixon
ProgramName = 'ortholog_outputter.py'
LastUpdated = '10/1/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Takes a set of orthologs and a list of fasta files and exports the 
sequences as fastas for alignment.

It reads a set of reciprocal orthologs, then reads the paired lists of protein and 
nucleotide fasta files to get ortholog sequences from. It then pulls 
the appropriate sequences for each ortholog from each fasta. 
The sequecnes for each chosen species are then output as fasta files
(one for the nucleotide one for the protein) named for the basal ortholog sequence
name, with sequence names in the fasta named for the species it came from.


'''

AdditionalProgramInfo = '''
Additional Program Information:
The protein and nucleotide fasta files need to have similar notation, each indicating the species
it came from in the same position in the name. For example:
Amillepora_Nuc.fasta and
Amillepora_Prot.fasta would work.

so would 
nuc_amil.fasta and
prot_amil.fasta

This wouldn't work
Amillepora_Nuc.fasta
amil_prot.fasta

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
parser.add_argument('-orthos', required = True, dest = 'orthos', help = 'The table of orthologs. Should be tab delimited and should have species names as the column headings that match the species names included in fasta file names')
parser.add_argument('-prot', required = True, dest = 'prot', nargs="+", help = 'A glob to the protein fastas, for example *PRO.fas')
parser.add_argument('-nucl', required = True, dest = 'nucl', nargs="+", help = 'A glob to the nucleotide fastas, for example *CDS.fas')
parser.add_argument('-spp', required = False, dest = 'spp', default = 0, help = 'Integer designating the pythonic position of the species identifier in the names of the fasta files. The default is 0, which would apply to a fasta file name like this: Amillepora_protein_seqs.fasta. Note that ortholog headings in the ortholog table must match the species identifiers in the fasta file names.')
parser.add_argument('-sep', required = False, dest = 'sep', default = '_', help = 'The delimiter used in the fasta file names. (Default is "_")')
args = parser.parse_args()

#Assign Arguments
Orthos = args.orthos
Prots = args.prot
Nucs = args.nucl
Spp = args.spp
Spp = int(Spp)
Sep = args.sep

sppList = []
for i in Nucs:
    sppList.append(i.split(Sep)[Spp])
print "\nNucleotide Fastas detected for the following species:"
for i in sppList:
    print i


def read_orthos(Orthos):
    """Function to read in the ortholog list and build 
    a set of nested dictionaries linking each species'
    ortholog to the basal ortholog name (assumed to be the first column of 
    the table)
    Output is like this:
    {BaseOrtholog1 : {spp1 : spp1ortho, spp2 : spp2ortho etc.}, BaseOrtholog2 : {spp1 : spp1ortho, etc.}
    
    It also makes the converse nested dictionary where baseOrthos are values and inidvidual species orthos are keys.
    {Species1 : {sppOrtho1 : base1, sppOrtho2 : base2}, Species2 : {sppOrtho1 : base 1, etc.} etc.}
    """
    with open(Orthos, 'r') as infile:
        base2xDict = {}
        x2baseDict = {}
        orthoList = []
        lineNumber = 0
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split('\t')
            if lineNumber == 1:
                sppList = line[1:]
                continue
            base = line[0]
            base2xDict[base] = {}
            orthoList.append(base)
            for i in range(len(sppList)):
                spp = sppList[i]
                ortho = line[1:][i]
                base2xDict[base][spp] = ortho
                try:
                    x2baseDict[spp][ortho] = base
                except KeyError:
                    x2baseDict[spp] = {}
    return orthoList,base2xDict, x2baseDict
        
def pull_seqs(orthoList, base2xDict, x2baseDict, fileSet):
    seqDict = {}
    for fasta in fileSet:
        #pull the species name from the fasta file name
        species = fasta.split(Sep)[Spp]
        #set up the nested dictionary for this sepcies
        seqDict[species] = {}
        #parse the fasta
        fasSeqs = SeqIO.parse(open(fasta), 'fasta')
        #iterate through the seqs
        for seq in fasSeqs:
            try:
                x2baseDict[species]
            except KeyError:
                continue
                exit("Species {} is missing for some reason".format(species))
            try:
                #try to pull the base ortholog connected with this sequence id
                base = x2baseDict[species][seq.id]
            except KeyError:
                #if you can't find a base ortholog for this one then it is not one of the reciprocal orthologs
                continue
            #if you found a connection with a base ortholog, then record the sequence in the sequence dictionary
            seqDict[species][base] = seq.seq.upper()
    return seqDict

def output(sppList, orthoList, seqDict, outSuffix):
    for base in orthoList:
        outFileName = base + outSuffix
        with open(outFileName, 'w') as out:
            counter = 0
            for species in sppList:
                counter += 1
                if counter == 1:
                    carrot = ">"
                else:
                    carrot = "\n>"
                try:
                    out.write("{}{}\n{}".format(carrot, species, seqDict[species][base]))
                except KeyError:
                    continue
            
                
        
    


orthoList, base2xDict, x2baseDict = read_orthos(Orthos)
# print "\n\nx2baseDict:"
# for i in x2baseDict.keys():
#     print i
# print base2xDict.keys()
# exit()
nucSeqDict = pull_seqs(orthoList, base2xDict, x2baseDict, Nucs)   
protSeqDict = pull_seqs(orthoList, base2xDict, x2baseDict, Prots)        
for i in nucSeqDict.keys():
    print nucSeqDict[i]
output(sppList, orthoList, nucSeqDict, "_nuc.fasta")
output(sppList, orthoList, protSeqDict, "_prot.fasta")

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))

