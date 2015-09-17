#MBD-seq_analysis1_distribution.R
#Groves Dixon
#1/16/15
#This script is the first of a set for analysis and figure generation
#from methylation, and substitution rate data from A.millepora and 
#other anthozoan species. This script generates a baseline set of 
#dataframes that are used in all subsequent scripts.  
#it is recommmended that you keep all input and output files in one
#working directory that you use for all the scripts.

#upload MBD-seq data and format it
setwd('/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory')
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")
mdat = read.table("methylDESeq_v6_8-3-15_clusted_dupRemoved.txt")
mdat$isogroup = rownames(mdat)
mdat = mdat[,c(2,8)]
colnames(mdat) = c('mbd.score', 'isogroup')
head(mdat)

#upload the iso2seq table which we'll need later
iso2seq = read.table("amil_seq2iso.tab", col.names = c('contig', 'isogroup'))

#upload CpG data, filter and calculate CpGoe
cdat = read.table("AmilleporaCpG.txt", header = T)
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
plot(mbd.score~cpgOE, data = cpg.dat, pch = 19, cex = 0.2)
plot(mbd.score~cpgOE, data = cpg.dat, pch = 19, cex = 0.1, ylim = c(-3, 3), xlim = c(0,1.5))

#do some stats on the correlation
lm1 = lm(mbd.score~cpgOE, data = cpg.dat)
summary(lm1)
abline(lm1, col = 'purple', lwd = 2)
library(Hmisc)
spear.cor = cor(cpg.dat$mbd.score, cpg.dat$cpgOE, method = "spearman")
spear.test = cor.test(cpg.dat$mbd.score, cpg.dat$cpgOE, method = "spearman")
spear.test


#look at the clustering of the MBD-scores
#this can take a little while, so it's set up 
#to run on a subset of the data. Comment in the 
#indicated line if you want to run it on all the MBD-scores

#set up a subset of the data
require(mclust)
y = mdat$mbd.score
x = sample(y, 1000)

#if you want to run it for real comment this line in and continue
x = mdat$mbd.score         #comment this back in to do actual figure

#use Bayesian Information Criterion to select the optimal mixture model/number of components
bic = mclustBIC(x, G = c(1:5))
summary(bic)
plot(bic, G = c(1:5))

#fit the optimal model
mod2 = Mclust(x, G = 3, modelNames = c('V'))

#make a plot showing the distribution of MBD-scores
#and the locations/densities of the mixture components
par(mfrow = c(2, 1))
par(mar = c(0, 4, 0, 2) + .1)
hist(mdat$mbd.score, breaks = 70, xlim = c(-5.3, 10), main = "", xlab = "\n", axes = F)
at = c(-5, 0, 5, 10)
axis(2, at = c(0, 1e3, 2e3), las = 1)
par(mar = c(5, 4, 0, 2) + .1)
plot(mod2, what = "classification", col = c('red', 'green', 'purple'), axes = F, main = "\n", xlab = "\n", xlim = c(-5.3, 10))
# abline(v = at)
axis(1, at = at)
title(ylab = 'Component')
title(xlab = "MBD-score")

#assemble the component assignments and merge with MBD dataset (only works for full dataset)
isogroup = mdat$isogroup
result = data.frame(isogroup, mod2$classification)
mdat2 = merge(mdat, result, by = 'isogroup')
colnames(mdat2) = c('isogroup', 'mbd.score', 'class')
head(mdat2)

#reassign the genes in component 3 to the strongly methylated or weakly methylated categories
mdat2[mdat2$class == 3 & mdat2$mbd.score > 0,]$class <- 2
mdat2[mdat2$class == 3 & mdat$mbd.score < 0, ]$class <- 1

#find the point where the two categories meet
strong = mdat2[mdat2$class == 2,]
weak = mdat2[mdat2$class == 1,]
max(weak$mbd.score)
min(strong$mbd.score)
separator = mean(c(max(weak$mbd.score), min(strong$mbd.score)))

#finalize
remove(strong)
remove(weak)
mdat <- mdat2

#--------- plot figure 1 ---------------
quartz()
#plotting variables
cex.lab = 1.25
xlim = c(0, 1.5)
cex.axis = 1.25
adj = .95
alpha = 0.05

#plot distribution
par(mar = c(4,3.5,1.5, 0) + 0.1)
par(mfrow = c(1,2))
x = hist(mdat$mbd.score, breaks = 70, xlim = c(-4, 5.05), main = "", axes = F, xlab = "\n", ylab = "\n", col = d, border = 'black', plot = F)
mids = x$mids
d <- mids > separator
d[d == TRUE] <- 'green'
d[d == FALSE] <- 'red'
hist(mdat$mbd.score, breaks = 70, xlim = c(-4, 5.05), main = "", axes = F, xlab = "\n", ylab = "\n", col = d, border = 'black', plot = T)
axis(1, cex.axis = cex.axis, at = c(-4, 0, 4))
axis(2, cex.axis = cex.axis, at = c(0, 1000, 2000))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Frequency', line = 2.5)

#plot correlation with CpGoe
par(mar = c(4,3.5,1.5,1) + 0.1)
r = plot.lm('mbd.score', 'cpgOE', cpg.dat, '\n', '\n', '\n', limits = T, c(-3.5, 4), ylim = c(.1, 1.3))
title(paste(paste('r =', r[1]), "****", sep = ""), line = -1, adj = adj, font.main = 3)
title(ylab = expression("CpG"["o/e"]), line = 2.25)
title(xlab = 'MBD-score', cex.lab = cex.lab, line = 2.5)
#----------------------------------------

#saving environment image for other scripts

#first get rid of some excess variables
remove(cg)
remove(cdat)
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
write.table(out, 'MBD-scores_for_GO_MWU.csv', quote = F, row.names = F, sep = ",")
#write out for KOGG
write.table(out, 'MBD-scores_for_KOGG_MWU.csv', quote = F, row.names = F, sep = ",")



