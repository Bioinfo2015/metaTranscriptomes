#!/usr/bin/env python
##build_paml_control.py
##written 6/26/14 by Groves Dixon
ProgramName = 'build_paml_control.py'
LastUpdated = '11/11/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Outputs a Paml control file for a given sequence file and a given tree file.
If necessary it collapses taxa in the tree file so that they match those given in 
the sequence file.

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
from Bio import Phylo
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input sequence file')
parser.add_argument('-t', required = True, dest = 'tree', help = 'The name of the tree file')
parser.add_argument('-spp', required = True, dest = 'sppList', help = 'A list of ALL species that could be included.')
parser.add_argument('-runMode', required = False, default = '-2', dest = 'runMode', help = 'The run mode you would like to use. The defulat is -2 (pairwise dN/dS calculation).')
args = parser.parse_args()




File = args.input
TreeFile = args.tree
SppFile = args.sppList
OutFileName = File.strip(".codon") + ".codeml"
ControlFileName = File.strip(".codon") + ".cnt"
TreeOutFileName = File.strip(".codon") + ".tree"
RunMode = args.runMode

def read_spp(SppFile):
    """read in the species list"""
    with open(SppFile, 'r') as infile:
        sppList = []
        for line in infile:
            sppList.append(line.strip("\n"))
    print "\nLooking for the following {} species in this file:".format(len(sppList))
    for i in sppList:
        print i
    return sppList

def read_codon_file(File, sppList):
    """read the codon file and extract the species that have
    sequences so that you can trim the tree file appropritately"""
    #set up lists to store the taxa that are present and absent in sequence file
    inSeqFileList = []
    absenteeList = []
    #read through file and find all the taxa 
    with open(File, 'r') as infile:
        for line in infile:
            line = line.strip("\n")
            if line in sppList:
                inSeqFileList.append(line)
    print "\nFound the following {} species in the sequence file:".format(len(inSeqFileList))
    for i in inSeqFileList:
        print i
    #build absentee list that will be removed from the phylogenetic tree to let codeml run for this sequence
    for i in sppList:
        if i not in inSeqFileList:
            absenteeList.append(i)
    print "\nThe following {} species were not found in the sequence file and will be removed from the tree:".format(len(absenteeList))
    for i in absenteeList:
        print i
    return absenteeList

def trim_tree(absenteeList, TreeFile):
    """Collapse away species from the phylogenetic tree that
    are not found in this sequence file. Output the tree file."""
    #parse the tree using Phylo
    tree = Phylo.read(TreeFile, 'newick')
    print "\nStarting with this Tree:\n"
    print(tree)
    #remove taxa that are not in the sequence file
    for taxon in absenteeList:
        tree.prune(taxon)
    print "\nExporting the this tree including only the species in this sequeunce file:\n"
    print(tree)
    #export the tree again
    Phylo.write(tree, TreeOutFileName, "newick")

#set up the base control file text for Codeml    
Base = '''
seqfile = {}
     treefile = {}
      outfile = {}

        noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = {}  * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

    cleandata = 0   * "I added on 07/07/2004" Mikita Suyama

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        model = 1
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = 0   * 0: simultaneous; 1: one branch at a time
'''.format(File, TreeOutFileName, OutFileName, RunMode)

sppList = read_spp(SppFile)
absenteeList = read_codon_file(File, sppList)
trim_tree(absenteeList, TreeFile)
#output the control file
with open(ControlFileName, 'w') as out:
    out.write(Base)


                
            
            

