#THIS SCRIPT USES DESeq TO ANALYZE MBDseq COUNTS
#1/15/15
#GROVES DIXON

########### IF YOU NEED TO DOWNLOAD DESeq2 #############
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")#upgrade biocLite if you need to
# biocLite("DESeq")
# biocLite("DESeq2")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
######################################################

#SET UP THE DATA TO RUN DESEQ
require('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/"
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/")

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF MATZ LAB RNA-SEQ PIPELINE. JUST APPLIED TO MBD DATA.
# counts=read.table('allCountsWithDupsCDS.txt',header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory. Older version
counts=read.table('mbd_raw_counts_8-1-15.txt', header=TRUE, row.names = 1) #these are mapped on the unclustered transcriptome without pcr dup removal
counts=read.table('clustered_premapping_counts.txt', header=TRUE, row.names = 1) #these are mapped on a clustered transcriptome without pcr dup removal
counts=read.table('get_old_dups_removed.txt', header=TRUE, row.names = 1) #these are mapped on a clustered transcriptome with pcr dup removal


head(counts)
# counts=read.table('allcounts_dupsRemoved__ubFisr_FIX.txt',header=TRUE,row.names=1) #NEW VERSION AS OF 5-12-15


# colnames(counts) = c("met1", "met2", "ub1", "ub2")
colnames(counts) = c("met1", "met2", "ub1", "ub2")
head(counts) 

#BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS
condition<-c("met", "met", "unbound", "unbound")
COLDATA <- data.frame(condition)

#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsHTSeq<-DESeqDataSetFromMatrix(counts, COLDATA, formula(~ condition))

#SET UP THE COLUMN DATA FOR THE DESeq TABLE
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("met","unbound"))
#note that putting 'untreated' first is so that direction of regulation in 'treated' group will be reflected in output

#RUN DESEQ
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)
tail(res)

##---------- LOOK AT THE RESULTS --------------------
#PLOT FOLD COVERAGE DIFFERENCES
quartz()
plotMA(dds,ylim=c(-4.5,9),main="DESeq2")

#LOOK AT THE DISTRIBUTION OF FOLD COVERAGES FOR ISOGROUPS
hist(res$log2FoldChange, breaks = 70)

##UPLOAD THE OLD CPGOE DATAFRAME AND MERGE WITH NEW DATA
cgm = read.table("/Users/grovesdixon/Documents/lab_files/projects/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores8-22-14__CDS.txt", header = T)
res$EST = rownames(res)
res2 = data.frame(res$EST, res$log2FoldChange)
colnames(res2) = c('EST', 'log2FoldChange')
head(res2)
dat = merge(res2, cgm, by = 'EST')
head(dat)

#PLOT CORRELATION BETWEEN FOLD CHANGE AND CPGOE
plot(log2FoldChange~cpgOE, data = dat, pch = 19, cex = 0.2)
lm1 = lm(log2FoldChange~cpgOE, data = dat)
abline(lm1, col = 'purple', lwd = 2)
summary(lm1)

plot(log2FoldChange~cpgOE, data = dat, pch = 19, cex = 0.1, ylim = c(-3, 3), xlim = c(0,1.5))
lm1 = lm(log2FoldChange~cpgOE, data = dat)
abline(lm1, col = 'purple', lwd = 2)



#LOOK AT CLUSTERING OF THE POINTS USING MCLUST
library(mclust)
data(diabetes)
class = diabetes$class
table(class)
X = diabetes[,-1]
head(X)
clPairs(X, class)
BIC = mclustBIC(X)
plot(BIC)
summary(BIC)
mod1 = Mclust(X)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")




head(dat)
y = data.frame(dat$log2FoldChange, dat$cpgOE)
colnames(y) = c('mbd', 'cpg')
z = y[1:500,]
head(z)
nrow(y)
head(y)
clPairs(y, one)
BIC = mclustBIC(z)
plot(BIC)
summary(BIC)
mod1 = Mclust(z, G = 3)
summary(mod1, parameters = TRUE)
mod1$G
mod1$BIC
quartz()
plot(mod1, what = "classification")
mod1$G
mod1$call
mod1$classification
z$class = mod1$classification
head(z)
plot(cpg~mbd, data = z, col = class)
#---------------------------------------------------------------

#PLOT CORRELATION FOR ONLY THE SIGNIFICANTLY METHYLATED AND NONMETHYLATED GENES
sig = dat[dat$pvalue < 0.05,]
plot(log2FoldChange~cpgOE, data = sig, pch = 19, cex = 0.2)
plot(log2FoldChange~CpG, data = sig)

#EXPORT THE DATASET
#USE MBD-seq_analysis.R TO ANALYZE AND BUILD FIGURES
write.table(res, "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v6_8-3-15_clusted_dupRemoved.txt", quote = F, sep = "\t")

