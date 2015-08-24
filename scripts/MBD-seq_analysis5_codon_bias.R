#MBD-seq_analysis5_codon_bias.R

#UPLOAD THE DATA
#mbd data
load("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")

#upload one of the codon bias datasets
cb = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/all_genes/Amillepora_RNA.out', header = T)


#format and merge with the mbd data
colnames(cb) = c('contig', 'cai', 'cbi', 'fop')
sub = merge(cb, iso2seq, by = 'contig')
cb = sub[,2:5]
head(cb)
head(mdat)
dat = merge(cb, mdat, by = 'isogroup')
head(dat)


#look at correlation between mbd.score and codon bias indices
#plot for using bias values based on all genes
QUANTILE = T
line.width = 1.5
quartz()
par(mfrow = c(2, 3))
xlim = c(-5, 6)
plot.lm('mbd.score', 'fop', dat, 'MBD-score', 'Frequency Optimal Codons', '\n', limits = T, xlim, c(.3, .7))
plot.lm('mbd.score', 'cbi', dat, 'MBD-score', 'Codon Bias Index', '\n', limits = T, xlim, c(-.25, .4))
plot.lm('mbd.score', 'cai', dat, 'MBD-score', 'Codon Adaptation Index', '\n', limits = T, xlim, c(.66, .85))

#plot error bar plots
do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, 'MBD-score', 'Frequency Optimal Codons', 'Fop', limits = F)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cbi', dat, 10, 'MBD-score', 'Codon Bias Index', 'CBI', limits = F)
plot.meth.bars(dat, 'mbd.score', 'cbi', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cai', dat, 10, 'MBD-score', 'Codon Adaptation Index', 'CAI', loess = F, limits = F)
plot.meth.bars(dat, 'mbd.score', 'cai', 0, c('green', 'red'))


#plot for using bias values based on all genes
QUANTILE = T
line.width = 2
cex.axis = 1.25
cex.lab = 1.25
quartz()
xlim = c(-5, 6)
plot.lm('mbd.score', 'fop', dat, '\n', '\n', '\n', limits = T, xlim, c(.3, .7))
plot.lm('mbd.score', 'cbi', dat, '\n', '\n', '\n', limits = T, xlim, c(-.25, .4))
plot.lm('mbd.score', 'cai', dat, '\n', '\n', '\n', limits = T, xlim, c(.66, .85))

#plot error bar plots
do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, '\n', '\n', '\n', limits = F)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cbi', dat, 10, '\n', '\n', '\n', limits = F)
plot.meth.bars(dat, 'mbd.score', 'cbi', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cai', dat, 10, '\n', '\n', '\n', loess = F, limits = F)
plot.meth.bars(dat, 'mbd.score', 'cai', 0, c('green', 'red'))


#plot same plots for using the top 5% highly expressed genes
QUANTILE = T
quartz()
par(mfrow = c(1, 1))
xlim = c(-5, 6)
line = -.5
adj = 1
#fop
r = plot.lm('mbd.score', 'fop', dat, '\n', '\n', '\n', limits = T, xlim, c(.3, .7))
title(paste(paste('r =', r[1]), "****", sep = ""), line = line, adj = adj, font.main = 1)
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Fop')
do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, '\n', '\n', '\n', limits = F)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Fop')

#cbi
plot.lm('mbd.score', 'cbi', dat, 'MBD-score', 'CBI', '\n', limits = T, xlim, c(-.25, .4))
do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, 'MBD-score', 'CBI', '\n', limits = F)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))

#cai
r = plot.lm('mbd.score', 'cai', dat, '\n', '\n', '\n', limits = T, xlim, c(.66, .85))
title(paste(paste('rho =', r[1]), "****", sep = ""), line = line, adj = adj, font.main = 1)
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'CAI')
do.window.plot(1/10, 'mbd.score', 'cai', dat, 10, '\n', '\n', '\n', loess = F, limits = F)
plot.meth.bars(dat, 'mbd.score', 'cai', 0, c('green', 'red'))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'CAI')


#plot the main figure
quartz()
par(mar = c(4,3.5,1.5,1) + 0.1)
par(mfrow = c(2, 2))
r = plot.lm('mbd.score', 'fop', dat, '\n', '\n', '\n', limits = T, xlim, c(.2, .8))
title(paste(paste('r =', r[1]), "****", sep = ""), line = line, adj = adj, font.main = 3)
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Fop', line = 2.5)
do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, '\n', '\n', '\n', limits = F)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Fop', line = 2.5)


r = plot.lm('mbd.score', 'cai', dat, '\n', '\n', '\n', limits = T, xlim, c(.66, .85))
title(paste(paste('r =', r[1]), "****", sep = ""), line = line, adj = adj, font.main = 3)
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'CAI', , line = 2.5)
do.window.plot(1/10, 'mbd.score', 'cai', dat, 10, '\n', '\n', '\n', loess = F, limits = F)
plot.meth.bars(dat, 'mbd.score', 'cai', 0, c('green', 'red'))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'CAI', , line = 2.5)



