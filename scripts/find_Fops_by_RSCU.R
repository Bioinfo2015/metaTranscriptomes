#find_Fops_by_RSCU
#find_Fops_by_RSCU
#This script generates a set of optimal codons based on RSCU in hi-expression and low-expression genes
#It outputs this list to generated Fop values for each coding sequence for A.millepora using calculate_frequency_optimal_codons.py

#set working directory
setwd("/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory")


#SET UP TOP 5 AND BOTTOM 5% OF GENES BASED ON EXPRESSION#SET UP TOP 5 AND BOTTOM 5% OF GENES BASED ON EXPRESSION
#keep this commented out if you already have them

#note that fewer genes in the low expression category have coding sequences
#extracted based on blastX. This is not surprising, since these genes will
#be harder to assemble in transcriptomes in general. To make up for this we
#grab the bottom 10%, which comes out to roughly the same number of contigs
#----------------------------------------------------
# sdat = read.table("mean_transcript_abundance.txt", header = T)
# sdat = sdat[order(sdat$rlog, decreasing = T),]
# head(sdat)
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
#----------------------------------------------------

#upload the top 5% and bottom 5% expressed gene lists
top5 = read.table('top5_hiExpressed_contigs.txt', header = T)
bottom5 = read.table('bottom5_lowExpressed_contigs.txt', header = T)
colnames(top5) = c('contig'); colnames(bottom5) = c('contig'); head(top5); head(bottom5)


#upload the RSCU data for all contigs
#upload the rscu data for all coding sequences from A.millepora
ru = read.table("Amillepora_all_RSCU.txt", header = T, na.strings = "NA")
#remove single use/start/stop codonds from this dataset
skip.codons = c('UAA', 'UAG', 'UGA', 'AUG', 'UGG')
ru = ru[, !(names(ru) %in% skip.codons)]
codons = colnames(ru)[2:length(colnames(ru))]
head(ru)

#gather means and perform t.tests between hi- and low-expressed genes for each codon
top = merge(ru, top5, by = 'contig')
bottom = merge(ru, bottom5, by = 'contig')
t.means = c()
b.means = c()
p.values = c()
for (c in codons){
	t = na.omit(top[,c])
	b = na.omit(bottom[,c])
	t.means = append(t.means, mean(t))
	b.means = append(b.means, mean(b))
	p = t.test(t, b, alternative = 'greater')$p.value
	p.values = append(p.values, p)
}
result = data.frame(codons, t.means, b.means, p.values)
result$adj.p = p.adjust(result$p.values, method = "BH")
head(result)

#assign optimal codons based on the p values from the t.tests
cut = 0.05
opt = result[result$adj.p < cut,]
opt1 = opt[order(opt$codons),]
nrow(opt1)

#check the overlap with the codonW optimal codons
opt.codonw = read.csv('~/Desktop/codonWopt.csv')
colnames(opt.codonw) = c('codons')
x = merge(opt1, opt.codonw, by = 'codons')
nrow(x)
print(paste('Percent overlap =', nrow(x)/nrow(opt)*100))


#output the optimal codon list for caluclating Fop with calculate_frequency_optimal_codons.py
#to get Fop values use this command:
#calculate_frequency_optimal_codons.py -i Amillepora_CDS.fas -optimal optimal_codons_by_rscu_hi-low_expressed.txt -o amil_Fop_by_rscu_hi-low.txt
head(opt1)
out = data.frame('aa', opt1$codons)
colnames(out) = c('AA', 'optimal')
write.table(out, "optimal_codons_by_rscu_hi-low_expressed.txt", row.names = F, quote = F)

#---------------------------------------------------------------------------------------------------------
#below are a couple other potential ways to assign optimal codons
#based on Spearman's rank correlations with expression that I decided not to use


#upload RSCU for all genes
ru = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/Amillepora_all_RSCU.txt", header = T, na.strings = "NA")

#upload expression data
sdat = read.table("mean_transcript_abundance.txt", header = T)
sdat = sdat[order(sdat$rlog, decreasing = T),]
head(sdat)
nrow(sdat)

#reformat the RSCU data
skip.codons = c('UAA', 'UAG', 'UGA', 'AUG', 'UGG')
#remove the stop and start codons from dataset
ru = ru[, !(names(ru) %in% skip.codons)]
codons = colnames(ru)[2:length(colnames(ru))]
ru = merge(ru, iso2seq, by = 'contig')
ru = merge(ru, sdat, by = 'isogroup')
head(ru)


#get Spearman's rho for each codon
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
result2 = data.frame(codons, rhos, ps)
result2$adj.p = p.adjust(result2$ps, method = "bonferroni")


#filter the dataset
pos = result2[result2$rho > 0,]
cut = 1e-10
opt2 = pos[pos$adj.p < cut,]

#check overlap with codonW optimal calls
nrow(opt2)
y = merge(opt2, opt.codonw, by = 'codons')
nrow(y)
print(paste('Percent overlap =', nrow(y)/nrow(opt2)*100))

#also try filtering directly on rho values rather than p values
r = pos$rhos
mean(r)
sd(r)
cut = mean(r) + sd(r)
plot(density(r))
abline(v = cut)
opt3 = result[pos$rho > cut,]
nrow(opt3)
z = merge(opt3, opt.codonw, by = 'codons')
nrow(z)
print(paste('Percent overlap =', nrow(z)/nrow(opt3)*100))




#or get overall usage data for concatenated forms of the hi and low expression datasets
hi.ru = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/top5_RSCU.txt', header = T, na.strings = "NA")
low.ru = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/codonW_analyses/get_RSCUs/bottom5_RSCU.txt', header = T, na.strings = "NA")
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
nrow(opt.dru)
opt.concensus = merge(opt.dru, opt.codonw, by = 'codons')
nrow(opt.concensus)
