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
directory<-"/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/noise/"
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/noise/")

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF MATZ LAB RNA-SEQ PIPELINE. JUST APPLIED TO MBD DATA.
counts1=read.table('allcounts.txt',header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory. #NOTE THIS TABLE WAS FIXED FOR LABELING
countsA = counts1[,1:3]
countsB = counts1[,7:9]
countsC = counts1[,13:15]
countsD = counts1[,19:21]
head(counts1)
counts2 = cbind(countsA, countsB, countsC, countsD)


#CRUDE STUFF ON RAW COUNTS
get.cv = function(vector){
  mn = mean(vector)
  std = sd(vector)
  cv = std/mn
  return(cv)
}
dat = countsA
dat$cv = apply(dat, 1, get.cv)
dat$mn = apply(dat, 1, mean)
dat = na.omit(dat)
head(dat)
dat = dat[with(dat, order(mn)),]
dat$log.cv = log(dat$cv)
dat$log.mn = log(dat$mn)
quartz()
plot(cv~log(mn, 10), data = dat, pch = 19, cex = .5)
loess_fit <- loess(cv~log(mn, 10), data = dat, span = .5, se = T)
lines(log(dat$mn, 10), predict(loess_fit),col="red",lwd=1)

#BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS
# condition<-c("a","a", "a", "a.s", 'a.s', 'a.s', 'b','b','b', 'b.s', 'b.s', 'b.s', 'c', 'c','c', 'c.s', 'c.s', 'c.s', 'd', 'd', 'd', 'd.s', 'd.s', 'd.s')
condition<-c('a','a','a', 'b', 'b','b','c','c','c','d','d','d')
COLDATA <- data.frame(condition)

#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsHTSeq<-DESeqDataSetFromMatrix(counts2, COLDATA, formula(~ condition))

#SET UP THE COLUMN DATA FOR TEH DESeq TABLE
# colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("a", "a.s", "b", "b.s", "c", "c.s", "d", "d.s"))
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('a','b','c', 'd'))
#note that putting 'untreated' first is so that direction of regulation in 'treated' group will be reflected in output

#RUN DESEQ
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)
tail(res)

plot(abs(log2FoldChange)~log(baseMean, 10), data = res)

#GET THE VARIANCE CONTROLLED EXPRESSION LEVEL
varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
ex = getVarianceStabilizedData(dds)
head(ex)
write.table(ex, "./variance_stabilized_expression_All_reps.txt")

ex$cv = apply(ex, 1, get.cv)
ex$mn = apply(ex[,1:3], 1, mean)
ex$cv = cv

v = merge(res, ex, by = intersect(names(y)))

plot(cv~log(mn, 10))
quartz()
plot(cv~log(mn, 10), data = d)
loess_fit <- loess(cv~log(mn, 10), data = d, span = .5, se = T)
lines(log(dat$mn, 10), predict(loess_fit),col="red",lwd=1)

#UPLOAD THE MBD-SEQ DATA
mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v1.txt")
mdat$EST = rownames(mdat)
head(mdat)


##---------- LOOK AT THE RESULTS --------------------
#PLOT FOLD COVERAGE DIFFERENCES
quartz()
plotMA(dds, ylim = c(-1, 1), main="DESeq2")

#PLOT RELATIONSHIP BETWEEN VARIATION AND MEAN EXPRESSION
plot(abs(log2FoldChange)~log(baseMean, 10), data = res)
z =res[with(res, order(baseMean)),]
plot(abs(log2FoldChange)~log(baseMean, 10), data = z)
loess_fit <- loess(abs(z$log2FoldChange)~log(z$baseMean, 10), span = 100, se = T)
lines(log(z$baseMean, 10), predict(loess_fit),col="red",lwd=1)
head(loess_fit)

#LOOK AT THE DISTRIBUTION OF FOLD COVERAGES FOR ISOGROUPS
hist(res$log2FoldChange, breaks = 70)

##UPLOAD THE OLD CPGOE DATAFRAME AND MERGE WITH NEW DATA
cgm = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores8-22-14__CDS.txt", header = T)
res$EST = rownames(res)
dat = merge(res, cgm, by = 'EST')

#PLOT CORRELATION BETWEEN FOLD CHANGE AND CPGOE
plot(log2FoldChange~cpgOE, data = dat, pch = 19, cex = 0.2)
plot(log2FoldChange~CpG, data = dat)