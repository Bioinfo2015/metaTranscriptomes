#MBD-seq_analysis6_positive_selection
#Groves Dixon
#8/28/15
#This script tests for relationships between MBD-score and positve selection
#as inferred from the branch-sites test for positve selection in PAML


#upload the data
setwd('/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory')
load("MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")
pdat = merge(read.table("branch_sites_LRT_results.txt", header = T), mdat, by = 'isogroup')
head(pdat)

#look at enrichment for positively selected genes
c = plot.counts(pdat, 'mbd.score', 'p.values', 1e-3, 1/20, 'MBD-score', 'Positive Selection', 0.5)
barplot(c$Ns)



#Get isogroup names for the contigs
head(iso2seq)
head(dnds)
sub = merge(dnds, iso2seq, by = 'contig')
head(sub)
dnds.dat = merge(sub, mdat, by = 'isogroup')[,1:7]
head(dnds.dat)

#make a preliminary plot for Adigitifera
require(plotrix)

