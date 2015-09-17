#LRT_for_branch_sites_models
#Groves Dixon
#8/28/15
#R script LRT_for_branch_site_models.R


#set working directory and 
setwd("/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory")
load("MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")

#read in the data from the null and alternative models for the branch-site test for positive selection
#see PAML manual for descriptions of these models
null = read.table("nullLikelihoods_branchSites.txt", header = T)
alt = read.table("altLikelihoods_branchSites.txt", header = T)

#these should show the contig, the number of parameters for the model, and it's log likelihood
head(null); head(alt)


#first run lrt between the models
contig = alt$contig
contigs.null = null$contig
la = alt$likelihood
lo = null$likelihood
p.values = lrt(la, lo, 1) #note degrees of freedom is difference in number of parameters between models (in this case always 1)

#record the results of the test in a new dataframe
result = data.frame(la, lo, p.values, contig, contigs.null)
head(result)

#merge the dataset with the contig-isogroup table 
dat = merge(result, iso2seq, by = 'contig')

#get adjusted p values for multiple tests
dat$adj.p = p.adjust(dat$p.values, method = 'BH')

#get an idea of how many genes are significant
cuts = c(0.05, 0.01, 0.001, 0.0001)
for (i in cuts){
	sub = dat[dat$adj.p < i,]
	unadjust = dat[dat$p.values < i,]
	print(paste(paste(i, nrow(sub)), nrow(unadjust)))
}


#write our the results
write.table(dat, "branch_sites_LRT_results.txt", quote = F, row.names = F)


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
write.table(out, 'branch_site_LRT_results_for_GOmwu.csv', quote = F, row.names = F, sep = ",")
#write out for KOGG
write.table(out, 'branch_site_LRT_results_for_GOmwu.csv', quote = F, row.names = F, sep = ",")