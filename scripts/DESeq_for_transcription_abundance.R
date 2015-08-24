#DESeq_for_transcription_abundance
#Analyze adult data from Dixon et al. 2015 to get transcript abundance estimates
#8/1/15
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
directory<-"/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/adultsVSjuvies/"
setwd(directory)

counts=read.table('adultVsLarval.txt', header=TRUE, row.names=1)
head(counts)

#SUBSET THE COUNTS DATA TO ONLY INCLUDE THE ADULTS THAT WERE NOT STRESS TREATED
sub = counts[,1:12]
counts <- sub
head(counts)

#EXPORT THE RAW READ COUNTS
raw = apply(counts, 1, mean)
out = data.frame(names(raw), raw)
head(out)
colnames(out) = c('isogroup', 'mn.count')
write.table(out, "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/raw_expression_counts.txt", quote = F, row.names = F)

#BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS
condition<-append(rep('oi', 6), rep('pcb', 6))
COLDATA <- data.frame(condition)

#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsHTSeq<-DESeqDataSetFromMatrix(counts, COLDATA, formula(~ condition))

#SET UP THE COLUMN DATA FOR THE DESeq TABLE
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("oi","pcb"))
#note that putting 'untreated' first is so that direction of regulation in 'treated' group will be reflected in output

#RUN DESEQ
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)
tail(res)

#PULL THE TOP 5% MOST HIGHLY EXPRESSED GENES
rld <- rlogTransformation(dds, blind=TRUE)
head(assay(rld))
x <- assay(rld)
mnExpression <- apply(x, 1, mean)
head(mnExpression)
five = round(length(mnExpression)/20)
mnExpression <- mnExpression[order(mnExpression)]
colnames(mnExpression) = c('isogroup', 'rlog')
head(mnExpression)
top.five = mnExpression[(length(mnExpression) - five):length(mnExpression)]
length(top.five)/length(mnExpression)
write.table(names(top.five), "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/high_expression/five_pct_highest_expressed_genes.txt", quote = F, row.names = F)
all.out = data.frame(names(mnExpression), mnExpression)
colnames(all.out) = c('isogroup', 'rlog')
write.table(all.out, "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/mean_transcript_abundance.txt", quote = F, row.names = F)

##---------- LOOK AT THE RESULTS --------------------
#PLOT FOLD COVERAGE DIFFERENCES
quartz()
plotMA(dds,ylim=c(-4.5,9),main="DESeq2")



#---------- OUTPUT RESULTS ----------------------
write.table(res, "./ageDiffExpression.txt", quote = F, sep = "\t")
