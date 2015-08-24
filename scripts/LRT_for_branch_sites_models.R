#R script LRT_for_branch_site_models.R
#read and analyze data from all five sites models in CODEML (model = 0, NSsites = 0, 1, 2, 7, 8)
#the different sites models are refered to by the value given for NSsites in the codeml control file
load("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")

#read in the data
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files")
null = read.table("nullLikelihoods.txt", header = T)
alt = read.table("altLikelihoods.txt", header = T)

#first run lrt between model 0 and model 1
#filter the dataset to only include genes that are somewhat significant
contig = alt$contig
contigs.null = null$contig
la = alt$likelihood
lo = null$likelihood
p.values = lrt(la, lo, 1) #note degrees of freedom is difference in number of parameters between models (in this case always 1)

#record the results of the test in a new dataframe
result = data.frame(la, lo, p.values, contig, contigs.null)
head(result)
dat = merge(result, iso2seq, by = 'contig')

#adjust the p values for multiple tests
dat$adj.p = p.adjust(dat$p.values, method = 'BH')

#get an idea of how many genes are significant
cuts = c(0.05, 0.01, 0.001, 0.0001)
for (i in cuts){
	sub = dat[dat$adj.p < i,]
	unadjust = dat[dat$p.values < i,]
	print(paste(paste(i, nrow(sub)), nrow(unadjust)))
}

#export the data for GO and KOGG enrichment tests using MWU-tests
#MISHA'S SCRIPTS USE ISOGROUP100 INSTEAD OF ISOGROUP=100, SO CHANGE THAT HERE
new.est = c()
for (i in dat$isogroup){
  y = paste(strsplit(i, '=')[[1]][1], strsplit(i, '=')[[1]][2], sep = '')
  new.est = append(new.est, y)
}
out = data.frame(new.est, -log(dat$p.values, 10))
colnames(out) = c('gene', 'logp')
head(out)
nrow(out)
#write out for GO
#use this file as input for 
write.table(out, '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/go/Acroporid_evolution_pvalues8-20-15.txt', quote = F, row.names = F, sep = ",")
#write out for KOGG
write.table(out, '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/kog/Acroporid_evolution_pvalues8-20-15.txt', quote = F, row.names = F, sep = ",")