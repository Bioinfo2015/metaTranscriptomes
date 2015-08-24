#TUTORIAL FOUND AT http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

########### IF YOU NEED TO DOWNLOAD DESeq2 #############
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")#upgrade biocLite if you need to
# biocLite("DESeq")
# biocLite("BiocGenerics")
# biocLite("DESeq2")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
######################################################

# require("BiocGenerics")
require('DESeq2')
#SAVE YOUR DIRECTORY NAME
directory<-"/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/pull_down/DESeq_folder/tutorial_dir/"

#GRAB YOUR TREATED FILES USING GREP ON 'TREATED' WHICH APPEARS IN ALL NAMES
sampleFiles <- grep("treated",list.files(directory),value=TRUE)

#BUILD A DATAFRAME ASSOCIATING EACH SAMPLE FILE WITH ITS TREATMENT CONDITION
sampleCondition<-c("treated","treated","treated","untreated","untreated","untreated")
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
head(sampleTable)

#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)

#SET UP THE COLUMN DATA
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))
#note that putting 'untreated' first is so that direction of regulation in 'treated' group will be reflected in output

#RUN IT
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)

#PLOT SOME RESULTS
quartz()
plotMA(dds,ylim=c(-2,2),main="DESeq2")

######################## REPEAT LOADING AS A FULL TABLE TO SEE THAT OUR WAY WILL WORK THE SAME #############
counts = read.table("/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/pull_down/DESeq_folder/tutorial_dir/singleTable.txt", header = T, row.names=1)
head(counts)

condition<-c("untreated", "untreated", "untreated","treated", "treated", "treated")
COLDATA <- data.frame(condition)

ddsHTSeq2<-DESeqDataSetFromMatrix(counts, COLDATA, formula(~ condition))
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))

#RUN IT
dds<-DESeq(ddsHTSeq2)
res<-results(dds)
res<-res[order(res$padj),]
head(res)

#PLOT SOME RESULTS
quartz()
plotMA(dds,ylim=c(-2,2),main="DESeq2")
############## PLOTS ARE IDENTICAL

