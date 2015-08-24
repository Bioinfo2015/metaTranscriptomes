#MBD-seq_analysis1_distribution.R
#Groves Dixon
#1/16/15
####################################################################################################
## THIS R SCRIPT ANALYZES MBD-seq DATA OUTPUT FROM DESeq (SEE DESeq_for_MBDseq.R)                 ##  
## GENES ARE LABELED AS ISOGROUPS CLUSTERED FROM THE A. MILLEPORA TRANSCRIPTOME (Moya et al. 2011)##
####################################################################################################

############# UPLOAD THE MBD FOLD COVERAGE RESULTS ############
#THE TABLE IS EXPORTED BY DESeq_ForMBDseq.R
setwd('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/')
# mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v1.txt")###previous dataset without duplicates removed
# mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v4_8-1-15.txt")#
# mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v5_8-3-15_clusted_dupIncluded.txt")#
mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v6_8-3-15_clusted_dupRemoved.txt")#
# mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v3_7-31-15.txt")##current dataset with duplicates removed. Mapped against transcriptome
mdat$isogroup = rownames(mdat)
head(mdat)

#upload the iso2seq table which we'll need later
iso2seq = read.table("/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/amillepora_transcriptome_july2014/amil_seq2iso.tab", col.names = c('contig', 'isogroup'))

###############################################################
############ LOOK AT CORRELATION WITH CPGOE ###################
###############################################################

#make new CpG data
cdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/cpgoe_data/AmilleporaCpG.txt", header = T)
head(cdat)
minlength = 100
maxlength = 20000 ###original value for gene bodies was 20000
minOE = .001  ##.001
maxOE = 2
cg = cdat
cg$cpgOE = (cg$CpG/(cg$C*cg$G))*(cg$length^2/(cg$length-1))##equation in Gavery and Roberts (oysters)
cg$gpcOE = (cg$GpC/(cg$C*cg$G))*(cg$length^2/(cg$length-1))
cg$tpgOE = (cg$TpG/(cg$T*cg$G))*(cg$length^2/(cg$length-1))
too.high = length(cg[cg$cpgOE > maxOE,]$EST)
too.low = length(cg[cg$cpgOE < minOE,]$EST)
too.short = length(cg[cg$length < minlength,]$EST)
print(paste("Too Short Transcript length", too.short))
print(paste("Too high CpG gene count =",too.high))
print(paste("Too low CpG gene count =",too.low))
cg = cg[cg$length > minlength,]
cg = cg[cg$length < maxlength,]
cg = cg[cg$cpgOE<maxOE,]
cg = cg[cg$cpgOE>minOE,]
cg = na.omit(cg)
head(cg)
colnames(cg)[1] = 'contig'

##UPLOAD THE OLD CPGOE DATAFRAME AND MERGE WITH NEW DATA
cgm = merge(cg, iso2seq, by = 'contig')
cpg.dat = merge(mdat, cgm, by = 'isogroup')

#PLOT CORRELATION BETWEEN FOLD CHANGE AND CPGOE
plot(log2FoldChange~cpgOE, data = cpg.dat, pch = 19, cex = 0.2)
plot(log2FoldChange~cpgOE, data = cpg.dat, pch = 19, cex = 0.1, ylim = c(-3, 3), xlim = c(0,1.5))

#DO STATS ON CORRELATION
lm1 = lm(log2FoldChange~cpgOE, data = cpg.dat)
summary(lm1)
abline(lm1, col = 'purple', lwd = 2)
library(Hmisc)
spear.cor = cor(cpg.dat$log2FoldChange, cpg.dat$cpgOE, method = "spearman")
spear.test = cor.test(cpg.dat$log2FoldChange, cpg.dat$cpgOE, method = "spearman")
spear.test

#############################################################
### LOOK AT CLUSTERING OF POINTS IN THE MBD VS CPGOE PLOT ###
#############################################################
#keep commented out unless want to plot figure 1
# library(cluster)
# x = data.frame(cpg.dat$cpgOE, cpg.dat$log2FoldChange)
# # #CLUSTER INTO TWO COMPONENTS
# colnames(x) = c('cpg', 'mbd')
# part <- pam(x, 2)
# id=part$clustering
# quartz()
# plot(mbd~cpg, data = x, col = id, ylim = c(-3, 3), xlim = c(0, 1.5), cex = .2, ylab = "log2 Fold Change", xlab = "CpGoe")
# # legend(1.05, 3.20, c('strong', 'weak'), fill = c('red', 'black'), title = "Methylation")
# legend(0, -2.3, c('strong', 'weak'), fill = c('red', 'black'), title = "Methylation")





#still in progress
## clustering using mclust
# library(mclust)
# head(dat)
# y = data.frame(cpg.dat$mbd.score, cpg.dat$cpgOE)
# colnames(y) = c('mbd', 'cpg')
# z = y[1:500,]
# head(z)
# nrow(y)
# head(y)
# # clPairs(y, one)
# BIC = mclustBIC(z)
# plot(BIC)
# summary(BIC)
# mod1 = Mclust(z, G = 3)
# summary(mod1, parameters = TRUE)
# mod1$G
# mod1$BIC
# quartz()
# plot(mod1, what = "classification")
# mod1$G
# mod1$call
# mod1$classification
# z$class = mod1$classification
# head(z)
# plot(cpg~mbd, data = z, col = class)











######################## PLOT FIGURE 1 ###################################
quartz()
hist(mdat$log2FoldChange, breaks = 70, xlim = c(5,-4), main = "", xlab = "Log2 Fold Difference")
abline(v = 0, lty = 2)
plot(mbd~cpg, data = x, col = c('darkgreen', 'black')[id], ylim = c(-3, 3), xlim = c(0, 2), cex = .3, ylab = "log2 Fold Change", xlab = "CpGoe", pch = 1)
# legend(1.31, -1.5, c('strong', 'weak'), fill = c('red', 'black'), title = "Methylation", bg = 'white')
legend(1.39, 3.215, c('strong', 'weak'), fill = c('black', 'darkgreen'), title = "Methylation", bg = 'white')
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
alpha.cols = add.alpha(c('black', 'darkgreen'), .1)
plot(mbd~cpg, data = x, ylim = c(-3, 3), xlim = c(0, 1.5), cex = .3, ylab = "log2 Fold Change", xlab = "CpGoe", pch = 19, col = alpha.cols[id])

###########################################################################
############ LOOK INTO POTENTIAL PULLDOWN-ONLY MEASURES ###################
###########################################################################

#UPLOAD THE COUNTS DATA
# setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/")
# dat = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/counts_2-10-15_dupsIncluded.txt', header = T)
# colnames(dat) = c("met1", "met2", "ub1", "ub2")#note this is correcting for the mislabeling before sequencing
# head(dat)
# 
# dat$EST = rownames(dat)

#################################################################
############# TEST FOR BIMODAL DISTRIBUTION #####################
#################################################################
#keep commented out 
# require(mclust)
# x = mdat$log2FoldChange
# # x = rnorm(n=500, m=1, sd=1) 
# # hist(x)
# clust = Mclust(x, G = c(1:5), modelNames = c("V"))
# bic = mclustBIC(x, G = c(1:9), modelNames = c("V"))
# plot(bic, G = c(1:9), modelNames = c("V"))

###########################################################
###### MODEL THE DISTRIBUTION WITH A MIXTURE MODEL ########
###########################################################

library(mixtools)
emmix = function(dat, col, means, sigma, lambda, k){
  x = dat[,col]
  cpgNull <- normalmixEM(x, mu = means, sigma = sigma, lambda = lambda, k=k, arbvar=T)
  summary(cpgNull)
  plot(cpgNull, which = 2, breaks = 50, density = F, xlab2 = "CpGo/e") 
  return(cpgNull)
}#function to find a best mixture model
mix.model = emmix(mdat, 'log2FoldChange', NULL, NULL, NULL, 2) #run the function for two components and no prior estimates of component parameters
model_df = function(lambda,mean,sigma){
  df = data.frame(cbind(lambda,mean,sigma))
  rownames(df) = paste("comp",rownames(df),sep="")
  return(df)
}##assigns the mixture model parameters to a table
mod.tab = model_df(mix.model$lambda, mix.model$mu, mix.model$sigma)

#############################################################################
########## PLOT THE DISTRIBUTION WITH THE MODEL COMPONENTS OVERLAID #########
#############################################################################

hist(mdat$log2FoldChange, breaks = 70, prob = T)
make_comp = function(model.name, comp.num){
  comp = function(x){
    model.name$lambda[comp.num]*dnorm(x, model.name$mean[comp.num], model.name$sigma[comp.num])
  }
  return(comp)
}#assigns the probability density function for each component for tracing
comp1 = make_comp(mod.tab, 1)
comp2 = make_comp(mod.tab, 2)
comp3 = make_comp(mod.tab, 3)
curve(comp1, add = TRUE, col='red', lwd=3)
curve(comp2, add = TRUE, col='green', lwd=3)
curve(comp3, add = TRUE, col='purple', lwd=3)
#OPTIONALL TRACE THE SUMMED MODEL 
pnormmix <- function(x,mixture){
  lambda <- mixture$lambda
  k <- length(mixture$lambda)
  pnorm.from.mix <- function(x,component){
    lambda[component]*dnorm(x,mean=mixture$mean[component],sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(rowSums(pnorms))
}
curve(pnormmix(x,mod.tab),add=T,col="purple",lty=1,lwd=2)


##save the data
mdat = mdat[,c(2,5:7)]
colnames(mdat) = c('mbd.score', 'pval', 'padj', 'isogroup')
cpg.dat = cpg.dat[,c(1,12:17)]
colnames(cpg.dat) = c('isogroup', 'GpC', 'TpG', 'length', 'cpgOE', 'gpcOE', 'tpgOE')
remove(cg)
remove(add.alpha)
remove(cdat)
remove(model_df)
save.image(file = "MBD-seq_Image.R")


#output the MBD-scores for GO and KOGG enrichment using MWU tests
#these output files are to be analyzed with GO_MWU_MBD-seq_enrichment.R and kog_MWU_MBD-seq_enrichment.R
new.est = c()
for (i in mdat$isogroup){
  y = paste(strsplit(i, '=')[[1]][1], strsplit(i, '=')[[1]][2], sep = '')
  new.est = append(new.est, y)
}
out = data.frame(new.est, mdat$mbd.score)
colnames(out) = c('gene', 'mbd.score')
head(out)
nrow(out)
#write out for GO
#use this file as input for 
write.table(out, '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/go/MBD-scores_for_GO_MWU.csv', quote = F, row.names = F, sep = ",")
#write out for KOGG
write.table(out, '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/kog/MBD-scores_for_KOGG_MWU.csv', quote = F, row.names = F, sep = ",")



