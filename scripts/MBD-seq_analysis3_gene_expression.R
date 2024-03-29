#MBD-seq_analysis3_gene_expression.R
#This script compares MBD-scores with patterns of gene expression
#Groves Dixon
#1/16/15

#upload R object and functions
setwd('/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory')
load('MBD-seq_Image.R')
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")

#upload the trascript abundance data 
sdat = read.table("mean_transcript_abundance.txt", header = T)
head(sdat)
head(cpg.dat)
edat = merge(mdat, sdat, by = 'isogroup')
e.cpg = merge(cgm, sdat, by = 'isogroup')

#plot relationship between mbd-score and transcript abundance

#set some plotting variables
require(plotrix)
QUANTILE = T
window = 1/25
cex.lab = 1.25
cex.axis = 1.25
line.width = 3
quartz()
par(mfrow = c(1, 2))
par(mar = c(4,3.5,1.5,.5) + 0.1)

#plot relationship with CpGoe
do.window.plot(window, 'cpgOE', 'rlog', e.cpg, 10, expression('\n'), '\n', '\n', F, .25, limits = T, xlim = c(1.25, .1), ylim = c(.4, 3))
axis(1, at = c(1.1, 0.7, 0.3), cex.axis = cex.axis)
axis(2, at = c(0.5, 1.75, 3), cex.axis = cex.axis)
title(xlab=expression("CpG"["o/e"]), ylab = 'Transcription', cex.lab = cex.lab, line = 2.5)

#plot for MBD-score
do.window.plot(window, 'mbd.score', 'rlog', edat, 10, '\n', '\n', '\n', F, .25, limits = F)
axis(1, at = c(-2, 0, 2), cex.axis = cex.axis)
axis(2, at = c(.5, 1.5, 2.5), cex.axis = cex.axis)
title(xlab= 'MBD-score', ylab = 'Transcription', cex.lab = cex.lab, line = 2.5)



#upload transplant gene expression data
transplantOrph = read.table("transplantOrphIn.txt", header = T)
transplantOrph$isogroup = rownames(transplantOrph) #so that it can be merged with cgm dataframe
head(transplantOrph)

#upload adults vs juveniles gene expression data
adat = read.table("ageDiffExpression.txt", header = T)
adat$isogroup = rownames(adat)
adat2 = data.frame(adat$isogroup, adat$padj, adat$log2FoldChange)
colnames(adat2) = c('isogroup', 'adj.p', 'log2.difference')

#merge the datasets with the mbd-seq data
mb = data.frame(mdat$isogroup, mdat$mbd.score)
colnames(mb) = c('isogroup', 'mbd.score')
ott = merge(mb, transplantOrph, by = 'isogroup')
avj = merge(mb, adat2, by = 'isogroup')

#take absolute values of log2 fold differences
ott$absDiff = abs(ott$log2.difference)
avj$absDiff = abs(avj$log2.difference)

#remove outliers
ott = ott[ott$absDiff < 7,]
oot = oot[oot$absDiff < 7,]
head(ott)
head(avj)


#----------- plot relationships between mbd-score and environmental variation in gene expression --------------
#the following code sets up a four panel figure relating MBD-score to the log2 fold differences and DEGs (p < 0.01) for environment and development

#set up some plotting variables
quartz()
par(mfrow = c(2,2))
par(mar = c(4,3.5,1.5,1) + 0.1)
title.line = -1.25
axis.line = 2.25
adj = .75
cex.lab = 1.25
cex.axis = 1.25
sep = separator
cut = 0.01

#development fold difference
z = cor.test(avj$mbd.score, avj$absDiff, method = "spearman")
rho = signif(z$estimate, digits = 2)
p = signif(z$p.value, digits = 2)
plot(absDiff~mbd.score, data = avj, xlim = c(-4, 6.2), col = alpha('grey', 0.2), ylim = c(0, 9), pch = 1, xlab = "\n", main = '\n', ylab = '\n', axes = F)
title(main = paste("r =", paste(rho, "****", sep = "")), adj = .75, line = title.line, adj = adj, font.main = 3)
title(xlab = "MBD-score", ylab = expression(paste("Fold Difference (log"['2'], ")", sep = "")), cex.lab = cex.lab, line = axis.line)
title(main = "Developmental")
axis(1, at = c(-3, 0, 3, 6), cex.axis = cex.axis)
axis(2, at = c(0, 4, 8), cex.axis = cex.axis)
lm1 = lm(absDiff~mbd.score, data = avj)
abline(lm1, col = 'red')
head(avj)

#development counts
plot.counts(avj, 'mbd.score', 'adj.p', 0.01, 1/25, '\n', '\n', .5)
title(xlab = "MBD-score", ylab = "DEGs", cex.lab = cex.lab, line = axis.line)
axis(1, at = c(4, 2, 0, -2), cex.axis = cex.axis)
axis(2, at = c(600, 750, 900), cex.axis = cex.axis)
x = do.meth.fisher(avj, sep, cut, 'mbd.score', 'adj.p', 'less')
p = x$p.value
odds = signif(x$estimate, digits = 2)
title(main = paste("OR =", paste(odds, "****", sep = "")), adj = .75, line = title.line, adj = adj, font.main = 3)
title(main = "Developmental")


#environmental fold difference
z = cor.test(ott$mbd.score, ott$absDiff, method = "spearman")
rho = signif(z$estimate, digits = 1)
p = signif(z$p.value, digits = 2)
res = paste(rho, p)
plot(absDiff~mbd.score, data = ott, xlim = c(-4, 6.2), col = alpha('grey', 0.2), ylim = c(0, 2), pch = 1, main = "\n", xlab = "\n", ylab = '\n', axes = F)
title(main = paste("r =", paste(rho, "****", sep = "")), adj = .75, line = title.line, adj = adj, font.main = 3)
title(xlab = "MBD-score", ylab = expression(paste("Fold Difference (log"['2'], ")", sep = "")), cex.lab = cex.lab, line = axis.line)
title(main = "Environmental")
axis(1, at = c(-3, 0, 3, 6), cex.axis = cex.axis)
axis(2, at = c(0, 1, 2), cex.axis = cex.axis)
lm1 = lm(absDiff~mbd.score, data = ott)
abline(lm1, col = 'red')


#environmental counts
plot.counts(ott, 'mbd.score', 'adj.p', 0.01, 1/25, '\n', '\n', .5)
title(xlab = "MBD-score", ylab = "DEGs", cex.lab = cex.lab, line = axis.line)
axis(1, at = c(4, 2, 0, -2), cex.axis = cex.axis)
axis(2, at = c(10, 20, 30, 40), cex.axis = cex.axis)
x = do.meth.fisher(ott, separator, cut, 'mbd.score', 'adj.p', 'less')
p = x$p.value
odds = signif(x$estimate, digits = 2)
title(main = paste("OR =", paste(odds, "****", sep = "")), adj = .75, line = title.line, adj = adj, font.main = 3)
title(main = "Environmental")

#---------------------------------------------------------------------

#make sanity plots for the above figure

#because the numbers of genes are so large, the stats can be misleading
#randomly reassign the MBD-score values and see what the correlations look like after that

#environmental sanity plots
s.ott = sanity(ott, 'mbd.score', 'absDiff')
z = cor.test(s.ott $mbd.score, s.ott $absDiff, method = "spearman")
rho = signif(z$estimate, digits = 2)
p = signif(z$p.value, digits = 2)
plot(absDiff~mbd.score, data = s.ott, xlim = c(-4, 6.2), col = alpha('grey', 0.2), ylim = c(0, 2), pch = 1, main = paste(rho, p), ylab = 'Development')
lm1 = lm(absDiff~mbd.score, data = s.ott)
abline(lm1, col = 'red')
s.ott = sanity(ott, 'mbd.score', 'adj.p')
plot.counts(s.ott, 'mbd.score', 'adj.p', 0.01, 1/25, 'MBD-score', 'Developmental DEGs', .5)
axis(1)
axis(2)
p = do.meth.fisher(s.ott, sep, cut, 'mbd.score', 'adj.p', 'less')
title(paste('p =', signif(p, digits = 2)))
box()

#developmental sanity plots
s.avj = sanity(avj, 'mbd.score', 'absDiff')
z = cor.test(s.avj$mbd.score, s.avj $absDiff, method = "spearman")
rho = signif(z$estimate, digits = 2)
p = signif(z$p.value, digits = 2)
plot(absDiff~mbd.score, data = s.avj, xlim = c(-4, 6.2), col = alpha('grey', 0.2), ylim = c(0, 9), pch = 1, main = paste(rho, p), ylab = 'Development')
lm1 = lm(absDiff~mbd.score, data = s.avj)
abline(lm1, col = 'red')
s.avj = sanity(avj, 'mbd.score', 'adj.p')
plot.counts(s.avj, 'mbd.score', 'adj.p', 0.01, 1/25, 'MBD-score', 'Developmental DEGs', .5)
axis(1)
axis(2)
p = do.meth.fisher(s.avj, sep, cut, 'mbd.score', 'adj.p', 'two.sided')
title(paste('p =', signif(p, digits = 2)))
box()

