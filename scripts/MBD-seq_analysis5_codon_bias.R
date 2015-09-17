#MBD-seq_analysis5_codon_bias.R
#This script analyzes the correlation between MBD-score and indices of codon bias
#It builds figures to show the correlations and tests the assumptions of the codonW analysis

#UPLOAD THE DATA
#mbd data
setwd("/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory")
load("MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")
library(plotrix)

#UPLOAD CODON BIAS DATASET
#from codonW
cb = read.table('codonW_results.txt', header = T)  #see codonW walkthrough for generating this file
head(cb)


#based on RSCU in hi- and low-expression genes
fop = read.table('amil_Fop_by_rscu_hi-low.txt', header = T)
head(fop)


#effective number of codons (also calculated by codonW)
nc = read.table("effectiveNumberCodons.txt", header = T)
head(nc)


#CAI based on exact calculations from top 5% highly expressed genes
cai = read.table("amil_exact_cai_values.txt", header = T)
head(cai)


#format and merge the codonW dataset with the mbd data
colnames(cb) = c('contig', 'cai', 'cbi', 'fop')
sub = merge(cb, iso2seq, by = 'contig')
cb = sub[,2:5]
head(cb)
head(mdat)
dat = merge(cb, mdat, by = 'isogroup')
head(dat)


#repeat for the fop dataset
sub = merge(fop, iso2seq, by = 'contig')
fop = merge(sub, mdat, by = 'isogroup')
head(fop)


#repeat for the Nc dataset
sub = merge(nc, iso2seq, by = 'contig')
nc = merge(sub, mdat, by = 'isogroup')
head(nc)


#repeat for the cai dataset
sub = merge(cai, iso2seq, by = 'contig')
cai = merge(sub, mdat, by = 'isogroup')
head(cai)

#look at correlation between mbd.score and all the different codon bias indices
QUANTILE = T
line.width = 1.5
quartz()
par(mfrow = c(2, 6))
xlim = c(-5, 6)
plot.lm('mbd.score', 'fop', dat, 'MBD-score', 'Frequency Optimal Codons', '\n', limits = T, xlim, c(.3, .7))
plot.lm('mbd.score', 'cbi', dat, 'MBD-score', 'Codon Bias Index', '\n', limits = T, xlim, c(-.25, .4))
plot.lm('mbd.score', 'cai', dat, 'MBD-score', 'CAI (codonW)', '\n', limits = T, xlim, c(.66, .85))
plot.lm('mbd.score', 'fop', fop, 'MBD-score', 'Frequency Optimal Codons', '\n', limits = T, xlim, c(.04,.21))
plot.lm('mbd.score', 'Nc', nc, 'MBD-score', 'Effective Number of Codons', '\n', limits = T, xlim, c(40, 61))
plot.lm('mbd.score', 'cai', cai, 'MBD-score', 'CAI exact', '\n', limits = T, xlim, c(0.6, 0.9))

#plot error bar plots
do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, 'MBD-score', 'Frequency Optimal Codons', 'Fop', limits = F)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cbi', dat, 10, 'MBD-score', 'Codon Bias Index', 'CBI', limits = F)
plot.meth.bars(dat, 'mbd.score', 'cbi', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cai', dat, 10, 'MBD-score', 'CAI', 'CAI(codonW)', loess = F, limits = F)
plot.meth.bars(dat, 'mbd.score', 'cai', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'fop', fop, 10, 'MBD-score', 'Frequency Optimal Codons', 'Fop', loess = F, limits = F)
plot.meth.bars(fop, 'mbd.score', 'fop', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'Nc', nc, 10, 'MBD-score', 'Effective Number of Codons', 'Nc', loess = F, limits = F)
plot.meth.bars(nc, 'mbd.score', 'Nc', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cai', cai, 10, 'MBD-score', 'CAI', 'CAI (exact)', loess = F, limits = F)
plot.meth.bars(cai, 'mbd.score', 'cai', 0, c('green', 'red'))

#----------- PLOT MAIN FIGURE ---------------------------
#plot the main figure based on CAI (codonW), Fop (RSCU hi-low), and Nc (codonW)
#first set up some plotting variables
require(plotrix)
quartz()
par(mar = c(3, 4.2, 1,1) + 0.1)
par(mfrow = c(3, 2))
cor.line = -1.6
cex.lab = 1.75
cex.axis = 1.5
lab.line = 2.75
cex.main = 1.5
adj = 1
line.width = 2
alpha = .075

#CAI
r = plot.lm('mbd.score', 'cai', cai, '\n', '\n', '\n', limits = T, xlim, c(.6, .9))
title(paste(paste('r =', r[1]), "****", sep = ""), line = cor.line, adj = adj, font.main = 3, cex.main = cex.main)
title(cex.lab = cex.lab, ylab = 'CAI', , line = lab.line)
wdat = do.window.plot(1/10, 'mbd.score', 'cai', cai, 10, '\n', '\n', '\n', loess = F, limits = F, xlim = c(-2.5, 3), ylim = c(0.7375, 0.781))
ydif = max(wdat$mn) - min(wdat$mn)
axis(1, cex.axis = cex.axis, cex.lab = cex.lab)
axis(2, at = c(0.73, 0.75, 0.77), cex.axis = cex.axis, cex.lab = cex.lab)
plot.meth.bars(cai, 'mbd.score', 'cai', 0, c('green', 'red'))
title(cex.lab = cex.lab, ylab = 'CAI', , line = lab.line)

#Fop
r = plot.lm('mbd.score', 'fop', dat, '\n', '\n', '\n', limits = T, xlim, ylim = c(.3, .7))
title(paste(paste('r =', r[1]), "****", sep = ""), line = cor.line, adj = adj, font.main = 3, cex.main = cex.main)
title(cex.lab = cex.lab, ylab = 'Fop', , line = 2.5)
wdat = do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, '\n', '\n', '\n', loess = F, limits = T, xlim = c(-2.5, 3), ylim = c(.45, .52))
ydif = max(wdat$mn) - min(wdat$mn)
axis(1, cex.axis = cex.axis, cex.lab = cex.lab)
axis(2, at = c(.45, .48, .51), cex.axis = cex.axis, cex.lab = cex.lab)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
title(cex.lab = cex.lab, ylab = 'Fop', , line = 2.5)

#Nc
par(mar = c(4.5, 4.2,.5,1) + 0.1)
# cor.line = -5
r = plot.lm('mbd.score', 'Nc', nc, '\n', '\n', '\n', limits = T, xlim, c(40, 61))
title(paste(paste('r = -0.30'), "****", sep = ""), line = cor.line, adj = adj, font.main = 3, cex.main = cex.main)
title(cex.lab = cex.lab, ylab = 'Nc', , line = 2.5)
title(xlab = 'MBD-score', line = 2.75, cex.lab = cex.lab)
wdat = do.window.plot(1/10, 'mbd.score', 'Nc', nc, 10, '\n', '\n', '\n', loess = F, limits = T, xlim = c(-2.5, 3), ylim = c(51.5,56))
ydif = max(wdat$mn) - min(wdat$mn)
axis(1, cex.axis = cex.axis, cex.lab = cex.lab)
axis(2, at = round(seq(from = min(wdat$mn), to = max(wdat$mn), by = ydif/2), digits = 0), cex.axis = cex.axis, cex.lab = cex.lab)
plot.meth.bars(nc, 'mbd.score', 'Nc', 0, c('green', 'red'))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Nc', , line = 2.75)
#--------------------------------------------------------------------------------------

#look at correlations with expression
sdat = read.table("mean_transcript_abundance.txt", header = T)
head(sdat)
c.e = merge(cai, sdat, by = 'isogroup')
f.e = merge(dat, sdat, by = 'isogroup')
n.e = merge(nc, sdat, by = 'isogroup')

quartz()
par(mfrow = c(3,1))
xlim = c(-4, 10)

#cai
r = plot.lm('rlog', 'cai', c.e, '\n', '\n', '\n', limits = T, xlim, c(.6, .9))
title(paste(paste(r[1]), "****", sep = ""), line = cor.line, adj = adj, font.main = 3, cex.main = cex.main)
title(cex.lab = cex.lab, ylab = 'CAI', , line = 2.5)
title(xlab = 'Expression', line = 2.75, cex.lab = cex.lab)
#fop
r = plot.lm('rlog', 'fop', f.e, '\n', '\n', '\n', limits = T, xlim, c(.3, .7))
title(paste(paste(r[1]), "****", sep = ""), line = cor.line, adj = adj, font.main = 3, cex.main = cex.main)
title(cex.lab = cex.lab, ylab = 'Fop', , line = 2.5)
title(xlab = 'Expression', line = 2.75, cex.lab = cex.lab)
#nc
r = plot.lm('rlog', 'Nc', n.e, '\n', '\n', '\n', limits = T, xlim, c(40, 61))
title(paste(paste(r[1]), "***", sep = ""), line = cor.line, adj = adj, font.main = 3, cex.main = cex.main)
title(cex.lab = cex.lab, ylab = 'Nc', , line = 2.5)
title(xlab = 'Expression', line = 2.75, cex.lab = cex.lab)




#---------- testing assumptions for automatic generation of Fop and CAI input files ----------
#upload expression data
sdat = read.table("mean_transcript_abundance.txt", header = T)
edat = merge(dat, sdat, by = 'isogroup')
edat = edat[edat$rlog != 0,]

#show correlation with expression

#fop
r = plot.lm('mbd.score', 'fop', dat, '\n', '\n', '\n', limits = F, xlim, c(.2, .8))
wdat = do.window.plot(1/10, 'rlog', 'fop', edat, 10, '\n', '\n', '\n', loess = F, limits = F)
axis(1); axis(2)

#cai
r = plot.lm('mbd.score', 'cai', dat, '\n', '\n', '\n', limits = F, xlim, c(.2, .8))
wdat = do.window.plot(1/10, 'rlog', 'cai', edat, 10, '\n', '\n', '\n', loess = F, limits = F)
axis(1); axis(2)

#make sure that top 5% biased genes are highly expressed
ax = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/all_genes/axis1.txt", header = T)
head(ax)
ax1 = merge(iso2seq, ax, by = 'contig')
ax = merge(sdat, ax1, by = 'isogroup')
head(ax)
ax = ax[with(ax, order(axis1)),]
head(ax)
p5 = round(0.05 * nrow(ax), 0)
top5 = na.omit(ax[1:p5,])
bottom5 = na.omit(ax[(nrow(ax) - p5) : nrow(ax),])
rest = na.omit(ax[p5+1 : nrow(ax) - p5,])

#set up the means and std.errors
t = mean(top5$rlog)
st = std.error(top5$rlog)
m = mean(rest$rlog)
sm = std.error(rest$rlog)
b = mean(bottom5$rlog)
sb = std.error(bottom5$rlog)
mns = c(t, m, b)
sterr = c(st, sm, sb)

#make dataframe for boxplots and barchart
top5$cat = 'top5'
rest$cat = 'middle'
bottom5$cat = 'bottom5'
bp = rbind(top5, rest, bottom5)

#plot the figures
boxplot(rlog~cat, data = bp)
barplot(mns, width = 1, space = .1, ylim = c(0, 2.7))
plotCI(c(0.6, 1.7, 2.8), mns, sterr, add = T, pch = 26, lwd = 2)
axis(1, at = c(0.6, 1.7, 2.8), labels = c("top 5%", "middle", "bottom 5%"), tick = F)
title(xlab = "Components of Principal Axis", ylab = "Mean Expression")

#test for correlation between principal axis and expression/effective number of codons

#expression
r = plot.lm('axis1', 'rlog', ax, '\n', '\n', '\n', limits = F, xlim, c(.2, .8))

#effective number of codons
nc1 = read.table("effectiveNumberCodons.txt", header = T)
head(nc1)
nc = merge(nc1, ax, by = 'contig')
head(nc)
r = plot.lm('axis1', 'Nc', nc, '\n', '\n', '\n', limits = F, xlim, c(.2, .8))

r = plot.lm('mbd.score', 'Nc', nc2, '\n', '\n', '\n', limits = T, xlim, c(.2, .8))

#look at enrichment for ribosomal genes
#upload the set of ribosomal isogroups
#these come from doing this: grep GO:0005840 amil_defog_iso2go.tab > ribosomal_isogroups.txt. (see CodonW_Instructions.txt for where to get these)
r = read.table("ribosomal_isogroups.txt", col.names = c('isogroup', 'go'))
r = read.table("glycolytic_isogroups.txt", col.names = c('isogroup', 'go'))


head(r)
head(top5)
nrow(top5)
head(ax)
rax = merge(ax, r, by = 'isogroup')
head(rax)

quartz()
par(mfrow = c(2,1))
xlim = c(-.6, 1)
hist(ax$axis1, xlim = xlim)
hist(rax$axis1, xlim = xlim)


#test for enrichment of ribosomal genes
p5 = round(0.05 * nrow(ax), 0)
top5 = na.omit(ax[1:p5,])
other = na.omit(ax[p5:nrow(ax), ])
nrow(top5)
nrow(others)
top.r = nrow(merge(top5, r, by = 'isogroup'))
other.r = nrow(merge(others, r, by = 'isogroup'))
top.not.r = nrow(top5) - top.r
other.not.r = nrow(other) - other.r

r.counts = c(top.r, other.r)
nr.counts = c(top.not.r, other.not.r)
m = as.table(rbind(r.counts, nr.counts))
colnames(m) = c('top5', 'other'); rownames(m) = c('r', 'not.r')
m
fisher.test(m)







abline(v = mean(ax$axis1))

head(dat)


xlim = c(0.6, 0.9)
hist(dat$cai, xlim = xlim)
rop = merge(dat, r, by = 'isogroup')
hist(rop$cai, xlim = xlim)


nrow(rop)

head(ax)
head(dat)


par(mfrow = c(3, 1))
hist(dat$cai, xlim = xlim)
meth = dat[dat$mbd.score > separator,]
not = dat[dat$mbd.score < separator,]
hist(meth$cai, xlim = xlim)
hist(not$cai, xlim = xlim)





