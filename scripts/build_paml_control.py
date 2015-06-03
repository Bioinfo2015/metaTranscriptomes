#!/usr/bin/env python
##build_paml_control.py
##written 11/11/14 by Groves Dixon
ProgramName = 'build_paml_control.py'
LastUpdated = '1/27/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Outputs a Paml control file for a given sequence file and a given tree file.

If necessary it collapses taxa in the tree file so that they match those given in 
the sequence file.

Now contains option to annotate the tree file  (adding "#1" to the branch branch leading to a particular clade 
for getting clade-specific dN/dS rates (see PAML_tree_tutorial_molevol2012 by Nicholas R. Meyerson for more on labeling trees for PAML)
'''

AdditionalProgramInfo = '''
Additional Program Information:
Default parameters for the model arguments are intended for running codeml in pairwise mode.

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
parser.add_argument('-model', required = False, default = '0', dest = 'model', help = 'The assignment for the model. Default is 0 (used for pair-wise comparisons).')
parser.add_argument('-NSsites', required = False, default = '0', dest = 'nssites', help = 'Set the values to use for NSsites. (Default = 0)')
parser.add_argument('-fix_omega', required = False, default = '0', dest = 'fixW', help = 'Set the value for fix_omega. (Default = 0)')
parser.add_argument('-omega', required = False, default = '0', dest = 'W', help = 'Set the fixed or initial value for omega. (Default = 0)')
parser.add_argument('-suffix', required = False, default = '', dest = 'suffix', help = 'A string to add to the output control file names so they can be differentiated from other control files for that gene. An example would be "AlternativeModel" or just "Alt"')
parser.add_argument('-treeNoteSuffix', required = False, default = '#1', dest = 'treeSuffix', help = 'A string that is added to terminal branches in the tree to designate the lineage for branch models in paml. (See paml manual and Tree_file_notes.txt for how the notations are added. Default value = "#1"')
parser.add_argument('-clade', required = False, default = 'none', dest = 'clade', help = 'A list of species (a subset of that given for -spp) that represent a clade within the tree. It is OK to include species that are missing from the tree, but the set of species must make up a monophyletic group within the tree. The given species will be used to assign the clade to measure dN/dS rates')
parser.add_argument('-inc', required = False, default = 'no', dest = 'inclusive', help = 'designate whether you want your "forground" dN/dS value to include terminal branches or just the branch leading to your clade of interest (Default = "no", type "yes" to be inclusive)')
args = parser.parse_args()


modelType = args.suffix
modelType = str(modelType)
File = args.input
TreeFile = args.tree
SppFile = args.sppList
if modelType != "":
    OutFileName = File.strip(".codon") + "_" + modelType + ".codeml"
    ControlFileName = File.strip(".codon") + "_" + modelType + ".cnt"
else:
    OutFileName = File.strip(".codon") + ".codeml"
    ControlFileName = File.strip(".codon") + ".cnt"
TreeOutFileName = File.strip(".codon") + ".tree"
CladeFile = args.clade
RunMode = args.runMode
Model = args.model
NSsites = args.nssites
FixOmega = args.fixW
Omega = args.W
SuffixInTree = args.treeSuffix
Inclusive = args.inclusive

if RunMode == "-2":
	TreeOutFileName = TreeFile
else:
	TreeOutFileName = File.strip(".codon") + ".tree"
import re


def read_spp_list(SppFile):
    """read in the species list"""
    with open(SppFile, 'r') as infile:
        sppList = []
        for line in infile:
            sppList.append(line.strip("\n"))
    return sppList

def read_codon_file(File, sppList):
    """read the codon file and extract the species that have
    sequences so that you can trim the tree file appropritately"""
    print "\nLooking for the following {} species in the sequence file:".format(len(sppList))
    for i in sppList:
        print i
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

def trim_tree(absenteeList, TreeFile, cladeList, Inclusive):
    """Collapse away species from the phylogenetic tree that
    are not found in this sequence file. Output the tree file."""
    print "\nReading the Tree..."
    #parse the tree using Phylo
    tree = Phylo.read(TreeFile, 'newick')
    terminals = tree.get_terminals()
    print "\nFound the following {} taxa in the tree:".format(len(terminals))
    for i in terminals:
        print i
    for taxon in absenteeList:
        tree.prune(taxon)
        if taxon in cladeList:
            cladeList.remove(taxon)
    print "\nAssigning clade in tree"
    print "\nThese species make up the clade:"
    for spp in cladeList:
        print spp
    #add "#1" to the branch leading to the clade of interest and all its terminal branches
    tree.common_ancestor(cladeList).name = "#1"
    #we may not want the terminal branches labeled. If so, strip the #1 from all terminal branch names
    if Inclusive == 'no':
        for leaf in tree.get_terminals():
            leaf.name = leaf.name.strip("#1")
    Phylo.write(tree, TreeOutFileName, "newick")
         
#set up the base control file text for Codeml
Base = '''
seqfile = {}
     treefile = {}
      outfile = {}

        noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = {}  * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

    cleandata = 0   * "I added on 07/07/2004" Mikita Suyama

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        model = {}
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = {}   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = {}   * 1: omega or omega_1 fixed, 0: estimate
        omega = {}   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = 0   * 0: simultaneous; 1: one branch at a time
'''.format(File, TreeOutFileName, OutFileName, RunMode, Model, NSsites, FixOmega, Omega)

sppList = read_spp_list(SppFile)
if CladeFile != "none":
    cladeList = read_spp_list(CladeFile)
absenteeList = read_codon_file(File, sppList)
#only bother trimming the tree if not runnin pairwise
if RunMode != '-2':
	trim_tree(absenteeList, TreeFile, cladeList, Inclusive)
#output the control file
with open(ControlFileName, 'w') as out:
    out.write(Base)


                
            
            

