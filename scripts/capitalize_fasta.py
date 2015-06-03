#!/usr/bin/env python
##capitalize_fasta.py
##written 9/16/14 by Groves Dixon
ProgramName = 'capitalize_fasta.py'
LastUpdated = '11/7/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:

Make all entries in a fasta capitalized


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
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-fa', required = True, dest = 'fa', help = 'The fasta file')
args = parser.parse_args()

#Assign Arguments
Fa = args.fa

     

fasSeqs = SeqIO.parse(open(Fa), 'fasta')
#iterate through the seqs
for seq in fasSeqs:
    print ">" + seq.id
    print seq.seq.upper()


#return time to run
Time = time.time() - Start_time
# print('\nTime took to run: {}'.format(Time))