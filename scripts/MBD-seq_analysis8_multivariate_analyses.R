#MBD-seq_analysis8_multivariate_analyses.R

#UPLOAD THE DATA
#mbd data
setwd("/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory")
load("MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")
library(plotrix)

#UPLOAD CODON BIAS DATASETS
#from codonW
cb = read.table('codonW_results.txt', header = T)  #see codonW walkthrough for generating this file
colnames(cb) = c('contig', 'CAI', 'CBI', 'Fop')
head(cb)
#based on RSCU in hi- and low-expression genes
fop = read.table('amil_Fop_by_rscu_hi-low.txt', header = T)
head(fop)
#effective number of codons (also calculated by codonW)
nc = read.table("effectiveNumberCodons.txt", header = T)
head(nc)

#UPLOAD EXPRESSION DATASET
sdat = read.table("mean_transcript_abundance.txt", header = T)
head(sdat)

#MERGE DATASETS
#fop
fop = merge(fop, iso2seq, by = 'contig')
fop = merge(fop, mdat, by = 'isogroup')
fop = merge(fop, sdat, by = 'isogroup')
head(fop)

#codon bias
cb = merge(cb, iso2seq, by = 'contig')
cb = merge(cb, mdat, by = 'isogroup')
cb = merge(cb, sdat, by = 'isogroup')
head(cb)

#effective number of codons
nc = merge(nc, iso2seq, by = 'contig')
nc = merge(nc, mdat, by = 'isogroup')
nc = merge(nc, sdat, by = 'isogroup')
head(nc)






#MULTIPLE REGRESSION
#fop
lm1 = lm(fop~mbd.score + rlog, data = fop)
lm1b = lm(fop~scale(mbd.score) + scale(rlog), data = fop)
summary(lm1b)


#codon bias (codonW)
lm2 = lm(CAI~mbd.score + rlog, data = cb)
lm2b = lm(CAI~scale(mbd.score) + scale(rlog), data = cb)
summary(lm2)
summary(lm2b)

#effective number of codons
lm3 = lm(Nc~mbd.score + rlog, data = nc)
lm3b = lm(Nc~scale(mbd.score) + scale(rlog), data = nc)
summary(lm3b)



#FOLLOWING ALONG WITH https://little-book-of-r-for-multivariate-analysis.readthedocs.org/en/latest/src/multivariateanalysis.html#multivariate-analysis

m = nc$mbd.score
var(m)
mean(m)
m2 = scale(m)
var(m2)
mean(m2)




