#MBD-seq_analysis5_codon_bias.R

#UPLOAD THE DATA
#mbd data
setwd("/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory")
load("MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")
library(plotrix)

#upload the codon bias datasets
cb = read.table('codonW_results.txt', header = T)  #based on codonW
rcb = read.table('amil_ribosome_Fop.txt', header = T, col.names = c('contig', 'total', 'optimal', 'nonoptimal', 'ambiguous', 'fop'))
head(cb)
head(rcb)

#format and merge the codonW dataset with the mbd data
colnames(cb) = c('contig', 'cai', 'cbi', 'fop')
sub = merge(cb, iso2seq, by = 'contig')
cb = sub[,2:5]
head(cb)
head(mdat)
dat = merge(cb, mdat, by = 'isogroup')
head(dat)

#repeat for the ribosomal dataset
sub = merge(rcb, iso2seq, by = 'contig')
rdat = merge(sub, mdat, by = 'isogroup')


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
plot.lm('mbd.score', 'fop', rdat, 'MBD-score', 'Frequency of Most Used', '\n', limits = T, xlim, c(.05,.28))

#plot error bar plots
do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, 'MBD-score', 'Frequency Optimal Codons', 'Fop', limits = F)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cbi', dat, 10, 'MBD-score', 'Codon Bias Index', 'CBI', limits = F)
plot.meth.bars(dat, 'mbd.score', 'cbi', 0, c('green', 'red'))
do.window.plot(1/10, 'mbd.score', 'cai', dat, 10, 'MBD-score', 'Codon Adaptation Index', 'CAI', loess = F, limits = F)
plot.meth.bars(dat, 'mbd.score', 'cai', 0, c('green', 'red'))


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


#plot the main figure
require(plotrix)
quartz()
par(mar = c(4,3.5,1.5,1) + 0.1)
par(mfrow = c(2, 2))
r = plot.lm('mbd.score', 'fop', dat, '\n', '\n', '\n', limits = T, xlim, c(.2, .8))
title(paste(paste('r =', r[1]), "****", sep = ""), line = line, adj = adj, font.main = 3)
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Fop', line = 2.5)
wdat = do.window.plot(1/10, 'mbd.score', 'fop', dat, 10, '\n', '\n', '\n', loess = F, span = 1, limits = F)
ydif = max(wdat$mn) - min(wdat$mn)
axis(1, cex.axis = cex.axis, cex.lab = cex.lab)
axis(2, at = signif(seq(from = min(wdat$mn), to = max(wdat$mn), by = ydif/2), digits = 3), cex.axis = cex.axis, cex.lab = cex.lab)
plot.meth.bars(dat, 'mbd.score', 'fop', 0, c('green', 'red'))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'Fop', line = 2.5)


r = plot.lm('mbd.score', 'cai', dat, '\n', '\n', '\n', limits = T, xlim, c(.66, .85))
title(paste(paste('r =', r[1]), "****", sep = ""), line = line, adj = adj, font.main = 3)
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'CAI', , line = 2.5)
wdat = do.window.plot(1/10, 'mbd.score', 'cai', dat, 10, '\n', '\n', '\n', loess = F, limits = F)
ydif = max(wdat$mn) - min(wdat$mn)
axis(1, cex.axis = cex.axis, cex.lab = cex.lab)
axis(2, at = signif(seq(from = min(wdat$mn), to = max(wdat$mn), by = ydif/2), digits = 3), cex.axis = cex.axis, cex.lab = cex.lab)
plot.meth.bars(dat, 'mbd.score', 'cai', 0, c('green', 'red'))
title(xlab = 'MBD-score', cex.lab = cex.lab, ylab = 'CAI', , line = 2.5)




#----------- look for codons with rscu that correlates with expression -----------
#upload the rscu data for all coding sequences from A.millepora
ru = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/Amillepora_all_RSCU.txt", header = T, na.strings = "NA")
sdat = read.table("mean_transcript_abundance.txt", header = T)
sdat = sdat[order(sdat$rlog, decreasing = T),]
head(sdat)
nrow(sdat)
skip.codons = c('UAA', 'UAG', 'UGA', 'AUG', 'UGG')
#remove the stop and start codons from dataset
ru = ru[, !(names(ru) %in% skip.codons)]
codons = colnames(ru)[2:length(colnames(ru))]
ru = merge(ru, iso2seq, by = 'contig')
ru = merge(ru, sdat, by = 'isogroup')
head(ru)




#set up top 5% and bottom 5% of genes based on expression
#note that fewer genes in the low expression category have coding sequences
#extracted based on blastX. This is not surprising, since these genes will
#be harder to assemble in transcriptomes in general. To make up for this we
#grab the bottom 10%, which comes out to roughly the same number of contigs
#keep this commented out if you already have them
#------------
# p5 = round(0.05 * nrow(sdat), 0)
# top5 = na.omit(sdat[1:p5,])
# bottom.bound = nrow(sdat) - p5*2
# bottom5 = na.omit(sdat[bottom.bound:nrow(sdat),])
# nrow(top5)
# nrow(bottom5)
# head(top5)
# head(bottom5)
# #write out the highest and lowest expressed genes
# t = merge(top5, iso2seq, by = 'isogroup')
# b = merge(bottom5, iso2seq, by = 'isogroup')
# write.table(t[,3], '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/top5_hiExpressed_contigs.txt', quote = F, row.names = F)
# write.table(b[,3], '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/bottom5_hiExpressed_contigs.txt', quote = F, row.names = F)
#--------------------

quartz()
par(mfrow = c(8, 8))
rhos = c()
ps = c()
for (c in codons){
	print(paste("Running Codon", c))
	u = ru[,c]
	e = ru$rlog
	t = na.omit(data.frame(u, e))
	# l = plot.lm('e', 'u', t, 'expression', 'usage', c, limits = F)
	z = cor.test(t[,'e'], t[,'u'], method = "spearman")
	rho = signif(z$estimate, digits = 2)
	p = z$p.value
	rhos = append(rhos, rho)
	ps = append(ps, p)
}
result = data.frame(codons, rhos, ps)
result$adjp = p.adjust(result$ps, method = "BH")
pos = result[result$rho > 0,]
x = pos$rhos
mean(result$rho)
sd(result$rho)
x = mean(result$rho) + sd(result$rho)
plot(density(result$rhos))
abline(v = x)

opt.rho = result[result$rho > x,]
nrow(opt.rho)

#-------------------- look at delta RSCU values -----------------------------------
hi.ru = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/top5_RSCU.txt', header = T, na.strings = "NA")
low.ru = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/bottom5_RSCU.txt', header = T, na.strings = "NA")
codons = names(hi.ru)[2:length(names(hi.ru))]
skip.codons = c('UAA', 'UAG', 'UGA', 'AUG', 'UGG')
hi = hi.ru[, !(names(hi.ru) %in% skip.codons)]
low = low.ru[, !(names(low.ru) %in% skip.codons)]
hi
low
codons = names(hi)[2:length(names(hi))]
hi = as.numeric(t(hi)[2:nrow(t(hi)), ])
low = as.numeric(t(low)[2:nrow(t(low)), ])
dru = data.frame(hi, low, row.names = codons)
dru$d = as.numeric(dru$hi) - as.numeric(dru$low)
plot(density(dru$d))
hist(dru$d)
mean(dru$d)
sd(dru$d)
cut = mean(dru$d) + sd(dru$d)
opt.dru = dru[dru$d > cut,]
opt.dru$codons = rownames(opt.dru)
opt.concensus = merge(opt.dru, opt.rho, by = 'codons')

opt.codonw = data.frame(c('UUU', 'UAU', 'UGU', 'UUA', 'UCA', 'CAU', 'CUA', 'CCA', 'CUG', 'AUU', 'AAU', 'AGU', 'AUA', 'ACA', 'AGA', 'AAG', 'AGG', 'GAU', 'AUA', 'GCA', 'GAA', 'GGA'))
colnames(opt.codonw) = c('codons')

c.dru = merge(opt.dru, opt.codonw, by = 'codons')
c.rho = merge(opt.rho, opt.codonw, by = 'codons')
c.dru = c.dru[order(c.dru$codons) ,]
c.rho = c.rho[order(c.rho $codons) ,]
nrow(c.dru)
nrow(c.rho)



#compare delta rscus between CG and GC codons
head(dru)
dru$codons = rownames(dru)
grep('CG', row.names(dru))
dru$cg = 0
dru$gc = 0
dru[grep('CG', dru$codons), "cg"] <- '1'
dru[grep('GC', dru$codons), "gc"] <- '1'
cg = na.omit(dru[dru$cg == 1,])
gc = na.omit(dru[dru$gc == 1,])
cg = cg[order(cg$codons),]
gc = gc[order(gc $codons),]
nrow(cg)
nrow(gc)
plot(density(cg$d))
lines(density(gc$d), col = 'red')
p = as.table(cbind(cg$d, gc$d))
names = append(rownames(cg), rownames(gc))
barplot(p, beside = T, names.arg = names, angle = 20)
abline(h = mean(cg$d))
abline(h = mean(gc$d))
counts <- table(mtcars$vs, mtcars$gear)


#barplot for codon adaptation index





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





