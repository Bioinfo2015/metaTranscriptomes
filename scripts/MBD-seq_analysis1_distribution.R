#MBD-seq_analysis1_distribution.R
#Groves Dixon
#1/16/15
#This script is the first of a set for analysis and figure generation
#from methylation, and substitution rate data from A.millepora and 
#other anthozoan species. This script generates a baseline set of 
#dataframes that are used in all subsequent scripts.  
#it is recommmended that you keep all input and output files in one
#working directory that you use for all the scripts.

#upload MBD-seq data
setwd('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/')
mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v6_8-3-15_clusted_dupRemoved.txt")#
mdat$isogroup = rownames(mdat)
head(mdat)

#upload the iso2seq table which we'll need later
iso2seq = read.table("/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/amillepora_transcriptome_july2014/amil_seq2iso.tab", col.names = c('contig', 'isogroup'))

#upload CpG data, filter and calculate CpGoe
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

#merge the CpGoe data with the MBD-seq data
cgm = merge(cg, iso2seq, by = 'contig')
cpg.dat = merge(mdat, cgm, by = 'isogroup')

#plot the correlation between mbd-score and CpGoe
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
require(mclust)
y = mdat$log2FoldChange
x = sample(y, 1000)
# x = mdat$log2FoldChange         #comment this back in to do actual figure

#
bic = mclustBIC(x, G = c(1:5))
summary(bic)
plot(bic, G = c(1:5))
isogroup = mdat$isogroup
mod2 = Mclust(x, G = 3, modelNames = c('V'))
par(mfrow = c(2, 1))
par(mar = c(0, 4, 0, 2) + .1)
hist(mdat$log2FoldChange, breaks = 70, xlim = c(-5.3, 10), main = "", xlab = "\n", axes = F)
at = c(-5, 0, 5, 10)
# abline(v = at)
# axis(1, at = at)
axis(2, at = c(0, 1e3, 2e3), las = 1)
par(mar = c(5, 4, 0, 2) + .1)
plot(mod2, what = "classification", col = c('red', 'green', 'purple'), axes = F, main = "\n", xlab = "\n", xlim = c(-5.3, 10))
# abline(v = at)
axis(1, at = at)
title(ylab = 'Component')
title(xlab = "MBD-score")
result = data.frame(isogroup, x, mod2$classification)
mdat2 = merge(mdat, result, by = 'isogroup')

head(result)
nrow(result)
head(cpg.dat)
r1 = data.frame(cpg.dat$isogroup, cpg.dat$cpgOE)
head(r1)
colnames(r1) = c('isogroup', 'cpgOE')
result2 = merge(result, r1, by = 'isogroup')
head(result2)
colnames(result2) = c('isogroup', 'mbd.score', 'class', 'cpgOE')
plot(mbd.score~cpgOE, data = result2, cex = .2)
plot(log2FoldChange~cpgOE, data = cpg.dat, pch = 19, cex = 0.2, col = result2$class)
plot(mbd.score ~cpgOE, data = result2, pch = 19, cex = 0.1, ylim = c(-3, 3), xlim = c(0,1.5), col = c(result2$class))












######################## PLOT FIGURE 1 ###################################
quartz()
hist(mdat$log2FoldChange, breaks = 70, xlim = c(-5, 10), main = "", xlab = "Log2 Fold Difference")
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
mix.model = emmix(mdat, 'log2FoldChange', NULL, NULL, NULL, 3) #run the function for two components and no prior estimates of component parameters
model_df = function(lambda,mean,sigma){
  df = data.frame(cbind(lambda,mean,sigma))
  rownames(df) = paste("comp",rownames(df),sep="")
  return(df)
}##assigns the mixture model parameters to a table
mod.tab = model_df(mix.model$lambda, mix.model$mu, mix.model$sigma)

#############################################################################
########## PLOT THE DISTRIBUTION WITH THE MODEL COMPONENTS OVERLAID #########
#############################################################################

hist(mdat$log2FoldChange, breaks = 70, prob = T, xlim = c(-5, 7))
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
curve(comp3, add = TRUE, col='yellow', lwd=3)
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

comp1



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



