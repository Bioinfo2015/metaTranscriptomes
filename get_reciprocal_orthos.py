#!/usr/bin/env python
##get_reciprocal_orthos.py.py
##written 9/16/14 by Groves Dixon
ProgramName = 'get_reciprocal_orthos.py'
LastUpdated = '10/21/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script pulls a set of reciprocal best hits based on 
a pair of reciprocal Blast outputs and the fasta files used as queries.
It has options to filter how the reciprocal blasts are selected.

'''

AdditionalProgramInfo = '''
Additional Program Information:

To use this script, first get two fasta files
that you think contain orthologs.
Blast each fasta against the other using out format 5 
so that biopython can parse the blast output file.
Now you have the two fasta files, and two blast output
files to pass as arguments to this script.

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
parser.add_argument('-br1', required = True, dest = 'br1', help = 'The first blast output file')
parser.add_argument('-br2', required = True, dest = 'br2', help = 'The second (reciprocal to first) blast output file.')
parser.add_argument('-fa1', required = True, dest = 'fa1', help = 'The fasta file that was the query for the first blast (and the database for the second blast)')
parser.add_argument('-fa2', required = True, dest = 'fa2', help = 'The fasta file that was the query for the second blast (and the database for the first blast)')
parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-id1', required = False, dest = 'id1', default = 0, type = int, help = 'The Position of the Identifying Information in the Sequence Names for fa1 (split by blank space). Default = 0')
parser.add_argument('-id2', required = False, dest = 'id2', default = 0, type = int, help = 'The Position of the Identifying Information in the Sequence Names for fa2 (split by blank space). Default = 0')
parser.add_argument('-e', required = False, dest = 'eval', default = 1e-20, help = 'Set an optional theshold to use for e-value. (Default is 1e-20)')
parser.add_argument('-hp', required = False, dest = 'hitPct', default = 80, type = float, help = 'Set an optional threshold for hit percentage. (Default 80)')
args = parser.parse_args()

#Assign Arguments
OutfileName = args.out
FA1 = args.fa1
FA2 = args.fa2
BR1 = args.br1
BR2 = args.br2
Ecut = args.eval
Ecut = float(Ecut)
ID1 = args.id1 ##the location of the identifying information in the seq names from fa1 when split by blank space
ID2 = args.id2 ##the location of the identifying information in seq names from fa2 when split by blank space
# ID1 = int(ID1)
# ID2 = int(ID2)
HitPctCut = args.hitPct / 100.0
NA = 'none'


# BR2 = args.br2
outfileName = args.out


def get_query_list(fa, Id):
    """Function to pull the full list of query sequences from a fasta file"""
    queryList = []
    for seq_record in SeqIO.parse(fa, "fasta"):
        queryList.append(seq_record.description.split()[Id])
    return queryList

def build_bestHit_dict(br, entryList, Id):
    """Builds a dictionary associating each entry from an fasta
    file with its best hit in the blast results file. Based on the 
    stringency of the blast, some entries will not have best hits.
    Their values are saved as 'none'. """
    print "Pulling the Best Hits for blast Results Files {}...".format(br)
    #first set up the dictionary
    hitDict = {}
    for i in entryList: hitDict[i] = 'none'
    #now read through the blast results to get the best hits
    result_handle = open(br) ##open the blast results file
    for blast_record in NCBIXML.parse(result_handle): ##start iterating through the blast records using the parser
        entry = blast_record.query
        entry = entry.split()[ID1]
        try:
            topHit = blast_record.alignments[0]
            hitName = topHit.title
            # print "entry {}     hitTitle {}".format(entry, hitName)
            hitName = hitName.split()[1].split()[Id]
            # print "entry {}     hitName {}".format(entry, hitName)
            #check that the hit percent is above cutoff, if it isn't then change hitName to 'Pct_Over'
            # if topHit.hsps[0].
            # print
            # print
            # print topHit.hsps[0].match
            # print topHit.hsps[0].expect
            # print topHit.length
            # print topHit.hsps[0].identities
            
            #check that the hit percentage is higher than the cutoff
            hitPct = float(topHit.hsps[0].identities)/float(topHit.length)
            # print hitPct
            if hitPct < HitPctCut:
                # print 'failed hit percent'
                hitName = 'Hit_Pct_to_low'
            #look at the evalue for the top hit
            if topHit.hsps[0].expect > Ecut:
                hitName = 'Eval_Over'
            # print "entry = {} : Hit = {} eval = {}".format(entry, hitName, topHit.hsps[0].expect)
        except IndexError:
            topHit = 'no hit'
        hitDict[entry] = hitName
    return hitDict
        
def pull_reciprocals(fa1List, fa2List, fa1Hits, fa2Hits):
    """Function to go through all the entries from
    each transcriptome and retain the reciprocal ortholog
    pairs. The data are stored in terms of the fa1Hits"""
    print "\nPulling Reciprocal Orthologs Based on Top Hits"
    ##I coded separate dictionaries at first, but then decided a results dict was best that stored all
    recipDict = {}
    missDict = {}
    resultDict = {}
    eFailList = []
    hpFailList = []
    noHitList = []
    for i in fa1List:
        bestHit = fa1Hits[i]
        #if the best hit was above the evalue threshold, store as an Efail and continue
        if bestHit == 'Eval_Over':
            eFailList.append(i) #note the evalue fail
            resultDict[i] = NA #record NA for the result
            continue
        if bestHit == 'Hit_Pct_to_low':
            hpFailList.append(i) #note the hit percentage fail
            resultDict[i] = NA #record NA for the result
            continue
        #check the reciprocal hit
        recipHit = fa2Hits[bestHit]
        if recipHit == i:
            recipDict[i] = bestHit #note the reciprocal hit
            resultDict[i] = bestHit #record result for output
        else:
            missDict[i] = bestHit #note that the reciprocal was different
            resultDict[i] = NA #record NA for the output
        Pass = recipHit == i
        # print "entry = {}  hit = {}  reciprocal = {}  pass = {}".format(i, bestHit, recipHit, Pass)
    ##print Summary Data for Pulling Reciprocals
    total = len(fa1List)
    recips = len(recipDict.keys())
    misses = len(missDict.keys())
    eFails = len(eFailList)
    hpFails = len(hpFailList)
    noHits = len(noHitList)
    sumTot = recips + misses + eFails + noHits + hpFails
    totOutput = len(resultDict.keys())
    print "\n{} Total Entries Found in {}".format(total, FA1)
    print "\nEvalue Cutoff was {}".format(Ecut)
    print "\nHit Percentage Cutoff was {}".format(HitPctCut)
    print "\n{} Entries had Reciprocal Best Hits with {}. Percent = {}%".format(recips, FA2, float(recips)/float(total))
    print "\n{} Entries had Hits, but Were not Reciprocal. Percent = {}%".format(misses, float(misses)/float(total))
    print "\n{} Entries hit at the original Evalue cutoff (used in the blast), but failed with the theshold given here".format(eFails)
    print "\n{} Entires hit but did not have high enough hit percentages and were not recorded".format(hpFails)
    print "\n{} Entries Never Hit".format(noHits)
    print "\nTotal Entires = {}; Total Accounted For Here = {}; Total to Output = {}".format(total, sumTot, totOutput)
    if total == sumTot == totOutput:
        print "\nNice, all Entries are Accounted For and Will be Output"
    else:
        print "\nHmmm, we lost some entires somewhere :( Recommend Checking somethign out"
    return resultDict, recipDict, missDict, eFailList, noHitList
    
def output_recips(resultDict, fa1List):
    """Outputs a table of results. For now they are 
    output in terms of the fa1 file entries"""
    print "\nOutputting Results..."
    with open(OutfileName, 'w') as out:
        header = "{}_Seqs\t{}".format(FA1, FA2)
        out.write(header)
        for i in fa1List:
            outString = "\n{}\t{}".format(i, resultDict[i])
            out.write(outString)

fa1List = get_query_list(FA1, ID1)
fa2List = get_query_list(FA2, ID2)
fa1Hits = build_bestHit_dict(BR1, fa1List, ID1)
fa2Hits = build_bestHit_dict(BR2, fa2List, ID2)
resultDict, recipDict, missDict, eFailList, noHitList = pull_reciprocals(fa1List, fa2List, fa1Hits, fa2Hits)
output_recips(resultDict, fa1List)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


