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
directory<-"/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/adultsVSjuvies/"
setwd(directory)

counts=read.table('adultVsLarval.txt', header=TRUE, row.names=1)
head(counts)

#BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS
condition<-append(rep('adult', 12), rep('juvenile', 30))
COLDATA <- data.frame(condition)

#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsHTSeq<-DESeqDataSetFromMatrix(counts, COLDATA, formula(~ condition))

#SET UP THE COLUMN DATA FOR THE DESeq TABLE
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("juvenile","adult"))
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


#---------- OUTPUT RESULTS ----------------------
write.table(res, "./ageDiffExpression.txt", quote = F, sep = "\t")
