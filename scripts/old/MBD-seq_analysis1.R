#MBD-seq_analysis.R
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
# mdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v2_5-12-15.txt")##previous dataset with duplicates removed. Mapped against transcriptome
xdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/methylDESeq_v2_6-1-15.txt")##current dataset with fix for met mislabel
xdat$EST = rownames(xdat)
mdat = data.frame(xdat$EST, xdat$log2FoldChange, xdat$pvalue, xdat$padj)
remove(xdat)
colnames(mdat) = c('EST', 'mbd.score', 'pavlue', 'padj')
head(mdat)

###############################################################
############ LOOK AT CORRELATION WITH CPGOE ###################
###############################################################

##UPLOAD THE OLD CPGOE DATAFRAME AND MERGE WITH NEW DATA
cgm = read.table("/Users/grovesdixon/Documents/lab_files/projects/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores8-22-14__CDS.txt", header = T)
cpg.dat = merge(mdat, cgm, by = 'EST')

#PLOT CORRELATION BETWEEN FOLD CHANGE AND CPGOE
plot(mbd.score~cpgOE, data = cpg.dat, pch = 19, cex = 0.2)
plot(mbd.score~cpgOE, data = cpg.dat, pch = 19, cex = 0.1, ylim = c(-3, 3), xlim = c(0,1.5))

#DO STATS ON CORRELATION
lm1 = lm(mbd.score~cpgOE, data = cpg.dat)
abline(lm1, col = 'purple', lwd = 2)
library(Hmisc)
pearson.cor = cor(cpg.dat$mbd.score, cpg.dat$cpgOE, method = "pearson")
pear.test = cor.test(cpg.dat$mbd.score, cpg.dat$cpgOE, method = "pearson")
pear.test

#############################################################
### LOOK AT CLUSTERING OF POINTS IN THE MBD VS CPGOE PLOT ###
#############################################################
#keep commented out unless want to plot figure 1
library(cluster)
x = data.frame(cpg.dat$cpgOE, cpg.dat$mbd.score)
# #CLUSTER INTO TWO COMPONENTS
colnames(x) = c('cpg', 'mbd')
part <- pam(x, 2)
id=part$clustering
quartz()
plot(mbd~cpg, data = x, col = id, ylim = c(-3, 3), xlim = c(0, 1.5), cex = .2, ylab = "log2 Fold Change", xlab = "CpGoe")
# legend(1.05, 3.20, c('strong', 'weak'), fill = c('red', 'black'), title = "Methylation")
legend(0, -2.3, c('strong', 'weak'), fill = c('red', 'black'), title = "Methylation")

######################## PLOT FIGURE 1 ###################################
quartz()
hist(mdat$mbd.score, breaks = 70, xlim = c(5,-4), main = "", xlab = "Log2 Fold Difference")
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
alpha.cols = add.alpha(c('black', 'blue'), .1)
plot(mbd~cpg, data = x, ylim = c(-3, 3), xlim = c(0, 1.5), cex = .3, ylab = "log2 Fold Change", xlab = "CpGoe", pch = 19, col = alpha.cols[id])

###########################################################################
############ LOOK INTO POTENTIAL PULLDOWN-ONLY MEASURES ###################
###########################################################################

#UPLOAD THE COUNTS DATA
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/")
dat = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/counts_2-10-15_dupsIncluded.txt', header = T)
colnames(dat) = c("met1", "met2", "ub1", "ub2")#note this is correcting for the mislabeling before sequencing
head(dat)

dat$EST = rownames(dat)

#################################################################
############# TEST FOR BIMODAL DISTRIBUTION #####################
#################################################################
#keep commented out 
# require(mclust)
# x = mdat$mbd.score
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
mix.model = emmix(mdat, 'mbd.score', NULL, NULL, NULL, 2) #run the function for two components and no prior estimates of component parameters
model_df = function(lambda,mean,sigma){
  df = data.frame(cbind(lambda,mean,sigma))
  rownames(df) = paste("comp",rownames(df),sep="")
  return(df)
}##assigns the mixture model parameters to a table
mod.tab = model_df(mix.model$lambda, mix.model$mu, mix.model$sigma)

#############################################################################
########## PLOT THE DISTRIBUTION WITH THE MODEL COMPONENTS OVERLAID #########
#############################################################################

hist(mdat$mbd.score, breaks = 70, prob = T)
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

############################################################################
################## PLOT BIOLOGICAL PROCESSES ###############################
############################################################################
par(mfrow = c(1,1))
require(plotrix)
#READ IN THE DATA
x = read.table('process2iso.txt', header = F)# see appendix 'setting up go term data' in Methylation_Pulldown_Analysis_Walkthrough.txt
#THE ORIGINAL TABLE HAS THE GO LISTINGS, WHICH ARE CUMBERSOME, BUILD NEW FRAME WITH JUST ISOGROUPS AND PROCESS NAME
pdat = data.frame(x$V1, x$V3)
colnames(pdat) = c("EST", "process")
head(pdat)
#GET THE PROCESS NAMES
ps = levels(pdat$process)
#SET UP EMPTY VECTORS TO STORE MEAN AND STANDARD ERROR IN FOR EACH BIOLOGICAL PROCESS
means = c()
st.errs = c()
Ns = c()
for (i in ps){
  sub = pdat[pdat$process == i,]
  comb = merge(sub, mdat, by = 'EST')
  comb$mbd.score = comb$mbd.score*-1#this is because reversing the X axis screws up plotCI()
  means = append(means, mean(comb$mbd.score))
  st.errs = append(st.errs, std.error(comb$mbd.score))
  Ns = append(Ns, nrow(comb))
}#populates the vectors with the means, Ns and standard errors
y = data.frame(ps, means, st.errs, Ns)
y = y[with(y, order(means)),]
par(mar=c(3,14,1,1))
plotCI(y$means, 1:length(y$means), uiw = y$st.errs, liw = y$st.errs, err="x", yaxt="n", pch=19, lwd=2, sfrac = .0075, axes = F, ylab = "")
axis(2, at = 1:nrow(y), labels = F, tick = T, xlim = c(.5, 1.3))
axis(1, at = c(-.5, 0, .5, 1), labels = c('0.5', '0.0', '-0.5', '-1.0'))
#Look at the GO term order
#Then type out nicely formated form
go.labels = c('DNA metabolism', 'Ribosome biogenesis', 'Translation', 'RNA metabolism', 'Transcription', 'Cell cycle and proliferation', "Stress response", "Transport", "Development", "Defense response", "Cell adhesion", "Signal transduction", "Response to stimulus", "Cell-cell signaling")
y.labels = paste(go.labels," (",y$Ns,")", sep = "")
mtext(y.labels, side = 2, line = .75, at = 1:nrow(y), outer=F,las=1, cex = 1.1)

#IN ADDITION TO PLOTTING THESE, WE CAN RUN MISHA'S GO ENRICHMENT TEST DIRECTLY ON THE LOG2FOLD CHANGES
#THIS REQUIRES THAT WE OUTPUT THE DATA IN A PARTICULAR FORMAT. THEN ANALYZE USING GO_MWU_MBD-seq_enrichment.R
out = data.frame(mdat$EST, mdat$mbd.score)
colnames(out) = c('gene', 'mbd.score')
head(out)
corrected.iso = c()
for (i in 1:length(out$gene)){
  
  corrected.iso = append(corrected.iso, paste(strsplit(as.character(out$gene[i]), "=")[[1]][1], strsplit(as.character(out$gene[i]), "=")[[1]][2], sep = ""))
}
out$gene = corrected.iso
write.table(out, "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/go/MBD_foldChangeForGO_V2_5-14-15.txt", sep = ',', quote = F, row.names = F)
#write one for KOGS to
write.table(out, "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/kog/MBD_foldChangeForKOG_V2_5-14-15.txt", sep = ',', quote = F, row.names = F)


#############################################################################
####### TEST HYPOTHESIS ABOUT TRANSITIONS FROM HI CPG TO LOW CPG ############
#############################################################################

head(cpg.dat)
head(pdat)
plot(mbd.score~cpgOE, data = cpg.dat)
abline(v = .65)
abline(h = 0)
abline(v = .8)
abline(h = 1)
cpg.cut = .65
mbd.cut = 0
new.sub = cpg.dat[cpg.dat$cpgOE > cpg.cut & cpg.dat$mbd.score > mbd.cut,]
mbd.methylated = cpg.dat[cpg.dat$mbd.score > 0,]
new.count = nrow(new.sub)
not.new.count = nrow(mbd.methylated) - new.count
ps

hi.cpg = cpg.dat[cpg.dat$cpgOE > cpg.cut,]
new = hi.cpg[hi.cpg$mbd.score > mbd.cut,]
never = hi.cpg[hi.cpg$mbd.score < mbd.cut,]
plot(mbd.score~cpgOE, data = cpg.dat)
points(mbd.score~cpgOE, data = new.meth, col = 'green')
points(mbd.score~cpgOE, data = never.meth, col = 'blue')
hi.cpg.count = nrow(hi.cpg)

#set up desired cutoffs
mbd.hi.cut = 0
mbd.low.cut = 0
cpg.hi.cut = .65
cpg.low.cut = .65
plot(mbd.score~cpgOE, data = cpg.dat, col = "grey75")
abline(v = c(cpg.hi.cut, cpg.low.cut))
abline(h = c(mbd.hi.cut, mbd.low.cut))


#set up four quadrant subsets
hi.cpg = cpg.dat[cpg.dat$cpgOE >= cpg.hi.cut,]
low.cpg = cpg.dat[cpg.dat$cpgOE < cpg.low.cut,]
hi.mbd = cpg.dat[cpg.dat$mbd.score >= mbd.hi.cut,]
low.mbd = cpg.dat[cpg.dat$mbd.score < mbd.low.cut,]
new = hi.cpg[hi.cpg$mbd.score > mbd.hi.cut,]
never = hi.cpg[hi.cpg$mbd.score < mbd.low.cut,]
old = cpg.dat[cpg.dat$mbd.score >= mbd.hi.cut & cpg.dat$cpgOE <= cpg.low.cut,]
miss = cpg.dat[cpg.dat$mbd.score < mbd.low.cut & cpg.dat$cpgOE < cpg.low.cut,]

#plot the quadrants to prove it worked
plot(mbd.score~cpgOE, data = cpg.dat, col = "grey75")
points(mbd.score~cpgOE, data = hi.mbd, col = 'red') #HI MBD
points(mbd.score~cpgOE, data = low.mbd, col = 'orange') #LOW MBD
points(mbd.score~cpgOE, data = hi.cpg, col = 'grey75') #HI CPG
points(mbd.score~cpgOE, data = low.cpg, col = 'grey75') #LOW CPG
points(mbd.score~cpgOE, data = new, col = 'green') #NEW = GREEN
points(mbd.score~cpgOE, data = never, col = 'blue') #NEVER = BLUE
points(mbd.score~cpgOE, data = old, col = 'black') #OLD = BLACK
points(mbd.score~cpgOE, data = miss, col = 'purple') #MISS = PURPLE

#set up the counts
never.counts = nrow(never)
new.counts = nrow(new)
old.counts = nrow(old)
miss.counts = nrow(miss)

total.counts = nrow(hi.cpg)
#run chi squares for enrichment of bio processes in the 'newly methylated' subset
i = "cellAdhesion"
results = c()
in.counts = c()
out.counts = c()
p.totals = c()
in.group = new
out.group = never
for (i in ps){
  #subset for the process
  p.sub = pdat[pdat$process == i,]
  #set up the in and outgroups
  p.in.sub = merge(p.sub, in.group, by = 'EST')
  p.out.sub = merge(p.sub, out.group, by = 'EST')
  #set up counts
  total.p.counts = nrow(p.sub)
  p.in.counts = nrow(p.in.sub)
  p.out.counts = nrow(meth.p.sub) - p.in
  not.p.in.counts = nrow(in.group) - p.in.counts
  not.p.out.counts = nrow(out.group) - p.out.counts
  process.counts = c(p.in.counts, p.out.counts)
  not.process.counts = c(not.p.in.counts, not.p.out.counts)
  M = as.table(rbind(process.counts, not.process.counts))
  colnames(M) = c("in", "out")
  rownames(M) = c("process", "not.process")
  print(i)
  print(M)
  z=fisher.test(M, alternative = "greater")
  results = append(results, z$p.value)
  in.counts = append(in.counts, p.in.counts)
  out.counts = append(out.counts, (total.p.counts - p.in.counts))
  p.totals = append(p.totals, total.p.counts)
}
finish = data.frame(ps, results, in.counts, out.counts, p.totals)
finish$ratio = finish$in.counts/finish$p.totals
finish = finish[with(finish, order(ratio)),]
finish$total = finish$in.counts + finish$out.counts

barplot(finish$ratio, names.arg = finish$ps)
plot(ratio~total, data = finish)
lm1 = lm(ratio~total, data = finish)
summary(lm1)
abline(lm1)
########################################################################################################
####################### PLOT RELATIONSHIP WITH EXPRESSION ##############################################
########################################################################################################
#UPLOAD THE EXPRESSION DATA
tdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/CpGoe_Project/Data-Analysis_files/data_from_transplant_experiment/VSD_lbRT_nodup_nomiso__reduced.csv", header =T)
sdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/noise/variance_stabilized_expression_All_reps.txt", header = T)
sdat$EST = rownames(sdat)
sdat$grandmeans = apply(sdat[,1:12], 1, mean)
sdat = merge(sdat, mdat, by = "EST")
head(sdat)
head(tdat)
#MERGE WITH THE MBD-SEQ DATASET
sdat$EST = rownames(sdat)
tdat = merge(tdat, mdat, by = "EST")
nrow(tdat)
tsig = tdat[tdat$pvalue < 0.05,]
##### SLIDING WINDOW FUNCTION #########
#----------- expression data -------------
require(plotrix)
#FUNCTION TO DIVIDE GENES INTO QUANTILES AND RETURN DATASET OF THE MEANS AND STANDARD ERRORS
window_plot = function(size, X, Y, dat, cut){
  windows = quantile(dat[,X], probs = seq(0, 1, by = size), na.rm = T)
  if (QUANTILE == F){
    width = max(dat[,X]) * size
    windows = seq(min(dat[,X]), max(dat[,X]), by = width)
  }
  print(windows)
  mn = c()
  x = c()
  N = c()
  sterr = c()
  cut.count = 0
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat[,X] >= left,]
    sub = sub[sub[,X] < right,]
    n = length(sub[,1])
    if (n < cut){
      cut.count = cut.count + 1
      next
    }
    mn = append(mn,mean(sub[,Y]))
    sterr = append(sterr,std.error(sub[,Y]))
    x = append(x, median(sub[,X]))
    N = append(N, n)
  }
  print(paste(cut.count, "Windows were cut Because they had below", cut, "genes"))
  plot_dat = data.frame(mn,x,sterr, N)
  return(plot_dat)
}##function to get the window data to plot th
Cutoff = 10
plot.meth.bars = function(dat, X, Y, cut, COLORS){
  m = dat[dat[,X] >= cut,]
  u = dat[dat[,X] < cut,]
  mns = c(mean(m[,Y]), mean(u[,Y]))
  meds = c(median(m[,X]), median(u[,X]))
  ses = c(std.error(m[,Y]), std.error(u[,Y]))
  print(meds)
  plotCI(meds, mns, uiw = ses, col = c('black', 'darkgreen'), add = T, pch = 19, lwd = 4, cex = .01)
}
############# PLOT EXPRESSION FIGURE ############################
quartz()
plot(grandmeans~mbd.score, data = sdat, col = 'grey', xlim = rev(range(sdat$mbd.score)))
window_data = window_plot(1/10, 'mbd.score', 'grandmeans', sdat, Cutoff) ## data for Carly's RNA seq
window_data = na.omit(window_data)##remove the windows that had no genes in them
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'black', col = 'black', lwd = 2, pch = ".", cex = .01, xlim = rev(range(window_data$x)))
#SWITCH BACK TO DEFAULT MAR AFTER PLOTTING BIOLOGICAL PROCESSES
par(mar = c(5, 4, 4, 2) + 0.1)
#PLOT THE POINTS
plot(mn~x,data=window_data, main="", pch = 1, cex = 1, axes = F, cex.lab = 1, ylim = NULL, xlim = rev(range(window_data$x)), ylab = "Variance Stabilized Expression")#use for tdat
axis(1)
axis(2)
#PLOT THE LOESS LINE
# loess_fit <- loess(mn ~ x, window_data, span = 1, se = T)
# lines(window_data$x, predict(loess_fit),col="red",lwd=1)
#PLOT THE ERROR BARS
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = F, scol = 'grey', col = 'grey', lwd = 2, pch = ".", cex = .01, xlim = rev(range(window_data$x)))
plot.meth.bars(sdat, 'mbd.score', 'grandmeans', 0, c('black', 'darkgreen'))



#plot cloud with sanity check overlaid
QUANTILE = F
plot(sanity~mbd.score, data = sdat, col = 'grey', xlim = rev(range(sdat$mbd.score)), xlab = "MBD-score", ylab = "Transcript Abundance")
window_data = window_plot(1/10, 'mbd.score', 'sanity', sdat, Cutoff) ## data for Carly's RNA seq
window_data = na.omit(window_data)##remove the windows that had no genes in them
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'black', col = 'black', lwd = 2, pch = ".", cex = .01, xlim = rev(range(window_data$x)))
x = sample(sdat$grandmeans, length(sdat$grandmeans))
sdat$sanity = x
window_data = window_plot(1/10, 'mbd.score', 'sanity', sdat, Cutoff) ## data for Carly's RNA seq
window_data = na.omit(window_data)##remove the windows that had no genes in them
# plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'blue', col = 'blue', lwd = 2, pch = ".", cex = .01, xlim = rev(range(window_data$x)))

#plot zooms
QUANTILE = T
window_data = window_plot(1/20, 'mbd.score', 'grandmeans', sdat, Cutoff) ## data for Carly's RNA seq
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = F, scol = 'black', col = 'black', lwd = 3, pch = 26, cex = .01, xlim = rev(range(window_data$x)), xlab = "mbd.score", ylab = "Transcript Abundance")
window_data = window_plot(1/20, 'mbd.score', 'sanity', sdat, Cutoff) ## data for Carly's RNA seq
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'grey75', col = 'grey75', lwd = 1, pch = 26, cex = .001, xlim = rev(range(window_data$x)))
window_data = window_plot(1/20, 'mbd.score', 'grandmeans', sdat, Cutoff) ## data for Carly's RNA seq
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'black', col = 'black', lwd = 3, pch = 26, cex = 1, xlim = rev(range(window_data$x)))
# plot.meth.bars(sdat, 'mbd.score', 'sanity', 0, c('black', 'darkgreen'))
# plot.meth.bars(sdat, 'mbd.score', 'sanity', 0, c('orange', 'orange'))

#################################################################



#PLOT SINGLE COLORED BOLD ERROR BARS FOR THE METHYLATED AND NONMETHYLATED GENES
weak = tdat[tdat$mbd.score < 0,]$grandmeans
strong = tdat[tdat$mbd.score > 0,]$grandmeans
mn.w = mean(weak)
st.w = std.error(weak)
mn.s = mean(strong)
st.s = std.error(strong)
mns = c(mn.w, mn.s)
sts = c(st.w, st.s)
Xs = c(-2, 1.5)
cols = c('red', 'green')
plotCI(Xs, mns, uiw = sts, add = T, lwd = 4, col = cols)


x = sample(sdat$mbd.score, length(sdat$mbd.score))
sdat$sanity = x
sdat = sdat[with(sdat, order(mbd.score)), ]
loess_fit <- loess(mbd.score~grandmeans, sdat, span = .2, se = T)
lines(sdat$mbd.score, predict(loess_fit),col="red",lwd=3)
plot(grandmeans~mbd.score, data = sdat, col = 'grey75')
plot(grandmeans~sanity, data = sdat, col = 'grey75')

y = y[with(y, order(means)),]

######################################
#THIS IS TOO HARD, 
# library(cluster)
# ?clusplot()
# ?pam
# x = data.frame(sdat$mbd.score, sdat$grandmeans)
# colnames(x) = c('mbd', 'exp')
# #CLUSTER INTO TWO COMPONENTS
# part <- pam(x, 3)
# cpg = part$data[,'cpg']
# mbd=part$data[,'mbd']
# y = data.frame(cpg, mbd)
# part2 <- pam(y, 2)
# id=part$clustering
# plot(mbd~cpg, col = id, ylim = c(-3, 3), xlim = c(0, 1.5), cex = .2)
# 
# plot(mbd.score~cpgOE, data = dat)
# clusplot(part, shade = F, color = F, add = F, labels = 0)
# ?clustplot.default
################################################









####################



########################################################################################################
####################### PLOT ENVIRONMENTAL VARIATION ###################################################
########################################################################################################
#UPLOAD THE TRANSPLANT GENE EXPRESSION DATA
setwd('/Users/grovesdixon/Documents/lab_files/projects/CpGoe_Project/Data-Analysis_files/cpg_OE')
originOrph = read.table("originOrphIn.txt", header = T)
transplantOrph = read.table("transplantOrphIn.txt", header = T)
originOrph$EST = rownames(originOrph) #so that it can be merged with cgm dataframe
transplantOrph$EST = rownames(transplantOrph) #so that it can be merged with cgm dataframe
head(originOrph)
head(transplantOrph)
#MERGE WITH THE PULLDOWN DATA AND SET UP NEW DATAFRAME
ott = merge(mdat, transplantOrph, by = 'EST')
oot = merge(mdat, originOrph, by = 'EST')
ott$absDiff = abs(ott$log2.difference)
oot$absDiff = abs(oot$log2.difference)
#THERE ARE A COUPLE OUTLIERS THAT MAKE THE ERROR BARS WONKY REMOVE THEM HERE
ott = ott[ott$absDiff < 7,]
oot = oot[oot$absDiff < 7,]
head(ott)
head(oot)
#SET THE MAR TO DEFAULT VALUES
par(mar = c(5, 4, 4, 2) + 0.1)

#FUNCTION TO PLOT COUNTS OF SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES
plot.counts = function(dat, X, Y, threshold, size){
  Ns = c()
  Xs = c()
  windows = quantile(dat[,X], probs = seq(0, 1, by = size))
  print(windows)
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat[,X] >= left,]
    sub = sub[sub[,X] < right,]
    sig.sub = sub[sub[,Y] < threshold,]
    Ns = append(Ns, nrow(sig.sub))
    Xs = append(Xs, mean(sub[,X]))
  }
  results = data.frame(Ns, Xs)
  print(results)
  plot(Ns~Xs, data = results, xlim = rev(range(results$Xs)), xlab = 'Strength of Methylation', ylab = 'Environmentally Responsive Genes', axes = F)
  loess_fit <- loess(Ns~Xs, results, span = 1, se = T)
  lines(results$Xs, predict(loess_fit),col="red",lwd=1)
  return(results)
}

#BUILD THE ENVIRONMENT PLOTS
quartz()
par(mfrow = c(2,1))
hist(mdat$mbd.score, breaks = 50, prob = F, xlim = c(4,-4), ylab = 'Genes', main = '', xlab = 'Strength of Methylation')

plot.counts(ott, 'mbd.score', 'adj.p', 0.01, 1/25)
axis(1, at = c(4, 2, 0, -2))
axis(2, at = c(10, 20, 30))
box()
window_data = window_plot(1/15, 'mbd.score', 'absDiff', ott, 10) ## data for Carly's RNA seq
plot(abs(log2.difference)~mbd.score, data = ott, xlim = rev(c(-3.5, 6.2)), col = 'grey', ylim = c(0, 2))
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'black', lwd = 1, pch = 19, cex = 0.01)
lm1 = lm(abs(log2.difference)~mbd.score, data = ott)
abline(lm1, col = 'red')
plot(abs(log2.difference)~mbd.score, data = ott, xlim = rev(c(-2,2)), ylim = c(.23,.37), col = 'grey')
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'black', lwd = 1, pch = 19, cex = 0.01)
#BULID THE POPULATION-WISE PLOTS
plot.counts(oot, 'mbd.score', 'adj.p', 0.01, 1/25)
window_data = window_plot(1/15, 'mbd.score', 'absDiff', oot, 10) ## data for Carly's RNA seq
plot(absDiff~mbd.score, data = oot, xlim = rev(c(-3.5, 6.2)), col = 'grey', ylim = c(0, 2))
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'black', lwd = 1, pch = 19, cex = 0.01)
lm1 = lm(abs(log2.difference)~mbd.score, data = oot)
abline(lm1, col = 'red')
plot(absDiff~mbd.score, data = oot, xlim = rev(c(-2,2)), ylim = c(.2,.37), col = 'grey')
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'black', lwd = 1, pch = 19, cex = 0.01)


#plot counts as barplot
counts = plot.counts(ott, 'mbd.score', 'adj.p', 0.01, 1/8)
head(counts)
barplot(rev(counts$Ns))
axis(1, at = c(4, 2, 0, -2))
axis(2, at = c(10, 20, 30))
box()

################################################################
############ CORRELATION WITH OPTIMAL CODON USAGE ##############
################################################################
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/")
odat = read.table('amil_Fop.txt', header = T)
head(odat)
iso2seq = read.table("/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/amillepora_transcriptome_july2014/amil_seq2iso.tab", col.names = c('EST', 'isogroup'))
head(iso2seq)
odat2 = merge(odat, iso2seq, by = 'EST')
head(odat2)
colnames(odat2) = c('contig', 'tot', 'optimal', 'suboptimal', 'ambiguous', 'Fop', 'EST')
head(mdat)
odat3 = merge(odat2, mdat, by = 'EST')
head(odat3)

quartz()
plot(Fop~mbd.score, data = odat3)
lm1 = lm(Fop~mbd.score, data = odat3)
summary(lm1)
abline(lm1, col = 'red')

par(mfrow = c(1,1))
window_data = window_plot(1/20, 'mbd.score', 'Fop', odat3, 10)
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = F, scol = 'black', lwd = 2, pch = 19, cex = 0.01, xlim = rev(range(window_data$x)), xlab = "MBD-score", ylab = "Frequency of Optimal Codons")
head(odat3)
met = odat3[odat3$mbd.score > 0,]
ub = odat3[odat3$mbd.score < 0,]
xs = c(median(met$mbd.score), median(ub$mbd.score))
mns = c(mean(met$Fop), mean(ub$Fop))
ses = c(std.error(met$Fop), std.error(ub$Fop))
plotCI(xs, mns, uiw = ses, add = T, scol = c('green', 'red'), col = c('green', 'red'), lwd = 2, pch = 19, cex = .5)

########### RELATE FOP TO EXPRESSION ####################
sdat = read.table("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/noise/variance_stabilized_expression_All_reps.txt", header = T)
sdat$EST = rownames(sdat)
sdat$grandmeans = apply(sdat[,1:12], 1, mean)
odat4 = merge(sdat, odat3, by = "EST")
head(odat4)

window_data = window_plot(1/20, 'grandmeans', 'Fop', odat4, 10)
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = F, scol = 'black', lwd = 2, pch = 19, cex = 0.01, xlim = rev(range(window_data$x)), xlab = "MBD-score", ylab = "Frequency of Optimal Codons")



########### SEE IF IT IS STILL TRUE WITHOUT RIBOSOMAL GENES #####################
#SET UP THE DATASET
rib = read.table("ribosomalGenes.txt") #grab the ribosomal gene set
colnames(rib) = c('EST')
odat3$ribosomal = odat3$EST %in% rib$EST
no.rib = odat3[odat3$ribosomal == FALSE,] #set up the dataset without them
rib2 = odat3[odat3$ribosomal == TRUE,]
#check that the numbers make sense
nrow(rib2)
nrow(no.rib)
nrow(odat3) - nrow(rib)

#plot the figure again
quartz()
par(mfrow = c(1,1))
window_data = window_plot(1/8, 'mbd.score', 'Fop', no.rib, 10)
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = F, scol = 'black', lwd = 2, pch = 19, cex = 0.01, xlim = rev(range(window_data$x)))
met = no.rib[no.rib$mbd.score > 0,]
ub = no.rib[no.rib$mbd.score < 0,]
xs = c(median(met$mbd.score), median(ub$mbd.score))
mns = c(mean(met$Fop), mean(ub$Fop))
ses = c(std.error(met$Fop), std.error(ub$Fop))
plotCI(xs, mns, uiw = ses, add = T, scol = c('green', 'red'), col = c('green', 'red'), lwd = 2, pch = 19)

head(window_data)


save.image(file = "MBD-seq_analysis1.R")



#################################################################################################################
################ DATA AND FUNCTIONS TO COMPARE PAIR-WISE SUBSTITUTION RATE VARIATION FOR A. MILLEPORA ###########
#################################################################################################################
#UPLOAD THE DNDS DATA FROM VARIOIUS OPTIONS
load("MBD-seq_analysis1.R")

dnds = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_Adigitifer.txt', header = T)#for digitifera
#UPLOAD THE ORTHOLOG DATA
orthos = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/Adigitifera_CDS_Ortholog_Table_e5_hp50_c2.txt', header = T)##best
#UPLOAD THE SEQ2ISO DATA FOR MILLEPORA
iso2seq = read.table("/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/amillepora_transcriptome_july2014/amil_seq2iso.tab", col.names = c('EST', 'isogroup'))
head(iso2seq)
head(orthos)
head(dnds)
head(mdat)
colnames(orthos)[1] = 'EST'

#FUNCTION TO MERGE THE DATAFRAMES UNTIL WE HAVE ALL THE IMPORTANT DATA ARRANGED TOGETHER
merge.dats = function(met.dat, dnds.dat, ortho.dat, iso2seq, mut.column, remove.outliers){
  sub = dnds.dat[dnds.dat$species == 'Amillepora',]
  sub2 = merge(sub, ortho.dat, by = 'EST')
  print(head(sub2))
  sub3 = data.frame(sub2$EST, sub2$Amillepora, sub2$dN, sub2$dS, sub2$dNdS)
  colnames(sub3) = c("digEST", "EST", "dN", "dS", "dnds")
  sub4 = merge(iso2seq, sub3, by = 'EST')
  colnames(sub4) = c('isotig', 'EST', 'digEST', 'dN', 'dS', 'dnds')
  head(sub4)
  mut.dat = merge(met.dat, sub4, by = "EST")
  result = mut.dat
  print(paste("Number of Genes before Outlier Removal = ", nrow(result)))
  if (remove.outliers == TRUE){
    result = mut.dat[mut.dat[,mut.column] < (mean(mut.dat[,mut.column]) + SDs*sd(mut.dat[,mut.column])), ]
    result$dnds = result$dN/result$dS
  }
  print(paste("Number of Genes after Outlier Removal = ", nrow(result)))
  return(result)
}
# SDs = .25
dnds.dat = merge.dats(mdat, dnds, orthos, iso2seq, 'dN', F)

#USE THE REFORMATTED DATA TO PLOT RELATIONSHIP BETWEEN PULLDOWN SCORE AND SUBSTITUTION RATES ##########
#FUNCTION TO PLOT A SLIDING WINDOW FOR GIVEN Y COLUMN AND X COLUMN
plot.mut = function(ddat, Ycol, Xcol, XLAB, YLAB, MAIN, rev, SEPERATOR){
  ddat = ddat[abs(ddat[,Ycol]) != Inf,] #get rid of instances where dNdS was Inf because dS was 0
  ddat = na.omit(ddat)
  window_data = window_plot(window.fraction, Xcol, Ycol, ddat, Cutoff)
  N = sum(window_data$N)
  TITLE = paste(MAIN, "N=",N)
  if (rev == TRUE){
    plot(mn~x,data=window_data, pch = 1, cex = 1, axes = F, cex.lab = 1, xlim = rev(range(window_data$x)), ylim = c(min(window_data$mn - window_data$sterr), max(window_data$mn + window_data$sterr)), ylab = XLAB, xlab = YLAB, main = MAIN)
  }
  if (rev == FALSE){
    plot(mn~x,data=window_data, pch = 1, cex = 1, axes = F, cex.lab = 1, xlim = range(window_data$x), ylim = c(min(window_data$mn - window_data$sterr), max(window_data$mn + window_data$sterr)), ylab = XLAB, xlab = YLAB, main = MAIN)
  }
  axis(1)
  axis(2)
  #PLOT THE ERROR BARS
  plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'grey', lwd = 0.5)
  #PLOT THE LOESS LINE
  loess_fit <- loess(mn ~ x, window_data, span = .9, se = T)
  lines(window_data$x, predict(loess_fit),col="red",lwd=1)
  #PLOT THE TWO FRACTIONS
  met = ddat[ddat[,Xcol] > SEPERATOR,]
  ub = ddat[ddat[,Xcol] < SEPERATOR,]
  mns = c(mean(met[,Ycol]), mean(ub[,Ycol]))
  print(mns)
  ses = c(std.error(met[,Ycol]), std.error(ub[,Ycol]))
  Xs = c(median(met[,Xcol]), median(ub[,Xcol]))
  if (rev == TRUE){
    COL = c('red', 'green')
  }
  if (rev == FALSE){
    COL = c('green', 'red')
  }
  plotCI(Xs, mns, uiw = ses, col = COL, add = T, lwd = 2, pch = 19)
  return(window_data)
}#ARGUMENTS: ddat = a dataframe with dn/ds data generated with 'merge.dats()'; Ycol = the column name for the data you want plotted on Y axis; Xcol = the column name for the data to use as X axis; XLAB = the desired X axis label; YLAB = desired Y label; MAIN = the title; rev = T/F variable for whether to reverse the X axis; SEPERATOR = the number that separates the methylated genes from nonmethylated.

#####################################################################################################################
######### PLOT RELATIONSHIP BETWEEN MBD-SEQ DATA AND PAIR-WISE SUBSTITUTION RATES WITH A.MILLEPORA ##################
####################################################################################################################
#----------------------- FUNCTION ----------------------------------------------
plot_list_of_species = function(speciesList){
  for (species in speciesList){
    print(species)
    dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
    print(head(dnds.dat))
    N = nrow(dnds.dat)
    title = paste(species, "N =", N)
    window.dat = plot.mut(dnds.dat, MUT.TYPE, METH.TYPE, MUT.TYPE, METH.TYPE, title, TRUE, 0)
  }
}
plot_list_of_speciesCpGoe = function(speciesList){
  for (species in speciesList){
    print(species)
    dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
    print(head(dnds.dat))
    dnds.dat = merge(cgm, dnds.dat, by = 'EST')
    N = nrow(dnds.dat)
    title = paste(species, "N =", N)
    window.dat = plot.mut(dnds.dat, MUT.TYPE, METH.TYPE, MUT.TYPE, METH.TYPE, title, FALSE, .65)
  }
}
#---------------------------------------------------------------------------------
#PLOT ACROPORIDS
METH.TYPE = 'mbd.score'
SDs = 2
MUT.TYPE = 'dN'
FILT = T
window.fraction = 1/10
quartz()
par(mfrow = c(1,4))
acro = c("Adigitifer", "Ahyacinthu", "Apalmata", "Atenuis")
METH.TYPE = 'mbd.score'
plot_list_of_species(acro)
METH.TYPE = 'cpgOE'
plot_list_of_speciesCpGoe(acro)

#PLOT OTHER CORALS
quartz()
par(mfrow = c(1,4))
non.acro = c('Pastreoide', 'Pcarnosus', 'Spistillat', 'Pdamicorni')
METH.TYPE = 'mbd.score'
plot_list_of_species(non.acro)
METH.TYPE = 'cpgOE'
plot_list_of_speciesCpGoe(non.acro)
#PLOT ANOMONES
quartz()
par(mfrow = c(1,4))
actin = c('Nvectensis', 'Apallida')
METH.TYPE = 'mbd.score'
plot_list_of_species(actin)
METH.TYPE = 'cpgOE'
plot_list_of_speciesCpGoe(actin)

########## RUN FOR dS ##############

SDs = 2
MUT.TYPE = 'dS'
METH.TYPE = 'mbd.score'
quartz()
par(mfrow = c(4,1))
#PLOT ACROPORIDS
plot_list_of_species(acro)
#PLOT OTHER CORALS
plot_list_of_species(non.acro)
#PLOT ANOMONES
plot_list_of_species(actin)

########## RUN FOR dNdS ##############
SDs = 2
MUT.TYPE = 'dnds'
plot_list_of_species(acro)
#PLOT OTHER CORALS
plot_list_of_species(non.acro)
#PLOT ANOMONES
plot_list_of_species(actin)

#-----------------------------------------------------------------------

########################################################################
########## PLOT A SUMMARY FIGURE WITH THE LOESS LINES STACKED ##########
########################################################################

#FUNCTION TO PLOT LOESS LINES
plot.mut.stackF = function(speciesList, mut.ylim, Ycol, Xcol, XLAB, YLAB, MAIN, rev, SEPERATOR, XLIM, COLORS){
#   par(xpd=F)
  i = 0
  for (species in speciesList){
    print(species)
    i = i + 1
    line.color = COLORS[i]
    sub.species = substr(species, 1, 10)
    dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", sub.species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
    print(head(dnds.dat))
    dnds.dat = dnds.dat[abs(dnds.dat[,MUT.TYPE]) < Inf,]##in some cases for the acroporids, dS was zero but dN wasn't, so filter Infs
    window_data = window_plot(1/8, Xcol, Ycol, dnds.dat, Cutoff)
    if (i < 2){
      plot(mn ~ x, type = "n", xlim = XLIM, ylim = mut.ylim, data = window_data, axes = F, xlab = "mbd.score", ylab = YLAB)
      axis(1)
      axis(2)
      loess_fit <- loess(mn ~ x, window_data, span = .9, se = T)
      lines(window_data$x, predict(loess_fit), lwd=2, col = line.color)
    }
    if (i >= 2){
      loess_fit <- loess(mn ~ x, window_data, span = .9, se = T)
      lines(window_data$x, predict(loess_fit), lwd=2.5, col = line.color)
    }
  }
#   par(xpd=T)
#   legendY = mut.ylim[1] + (mut.ylim[2] - mut.ylim[1])*4/5
#   legend(-2.25, legendY, rev(speciesList), fill = rev(COLORS))
#   par(xpd=F)
}

sub.species = 'Adigitifer'
dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", sub.species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
print(head(dnds.dat))
dnds.dat = dnds.dat[dnds.dat[,MUT.TYPE] != Inf,]


######### PLOT FOR ACROPORID COMPARISONS ###############
#dN
quartz()
par(mfrow = c(1, 3))
tot = sum(c(length(acro), length(non.acro), length(actin)))
COLOR.SET = rainbow(tot)
COLS = COLOR.SET[1:4]
MUT.TYPE = 'dN'
MUT.YLIM = c(0.006, 0.025)
YLAB = "dN"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dS'
MUT.YLIM = c(.025, .125)
YLAB = "dS"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dNdS
MUT.TYPE = 'dnds'
MUT.YLIM = NULL
YLAB = "dN/dS"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

############### PLOT FOR OTHER CORALS ###############
#dN
###PLOT dN FOR NON ACROPORID SPECIES
quartz()
par(mfrow = c(1, 3))
COLS = COLOR.SET[5:8]
MUT.TYPE = 'dN'
MUT.YLIM = NULL
YLAB = "dN"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dS'
MUT.YLIM = c(1.8, 3.2)
YLAB = "dS"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dN/dS
MUT.TYPE = 'dnds'
MUT.YLIM = c(0.048, 0.09)
YLAB = "dN/dS"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

########## PLOT dS FOR ACTINIARIA #################
#dN
COLS = COLOR.SET[9:10]
quartz()
par(mar=c(5, 4, 0, 8) + 0.1)
par(xpd=F)
MUT.TYPE = 'dN'
MUT.YLIM = c(.185, .235)
YLAB = "dN"
plot.mut.stackF(actin, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
quartz()
par(mar=c(5, 4, 0, 8) + 0.1)
par(xpd=F)
MUT.TYPE = 'dS'
MUT.YLIM = c(24.5, 30)
YLAB = "dS"
plot.mut.stackF(actin, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

#PLOT THE LOESS LINES FOR PAIR-WISE MUTATION RATES FOR THE ACROPORID SPECIES ########

#########Plotting Parameters 
quartz()
par(mfrow=c(1,1))
par(mar=c(5, 4, 4, 8) + 0.1)
par(xpd=F)
XLAB = "mbd.score"
acroporids <- c("Adigitifera", "Ahyacinthus", "Apalmata", "Atenuis")
FILT = T

###dN specific plotting parameters
MUT.TYPE = 'dN'
MUT.YLIM = c(0.004, 0.016)
YLAB = "dN"
plot.mut.stackF(acroporids, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2))

###dS specific plotting parameters
MUT.TYPE = 'dS'
MUT.YLIM = c(.025, .125)
YLAB = "dS"
plot.mut.stackF(acroporids, MUT.YLIM, MUT.TYPE, 'mbd.score', XLAB, YLAB, species, TRUE, 0, c(2, -2))

#################### EXPORT DNDS FOR KOGG AND GO #########################
#here I want to run Kogg and GO on the dNdS values from each pairwise comparison
#The idea is that by looking at the types of genes that pop up in one type of comparison vs another could
#show what types of genes were under selection in the Acroporid lineage
#the following loop outputs a bunch of dnds tables for each species pair to the kog folder
#use these as input for the KOGG MWU script to look at the types of genes that evolve rapidly between pairs
spp.list = c("Adigitifer", "Ahyacinthu", "Apalmata", "Atenuis", "Pastreoide", "Pcarnosus", "Pdamicorni", "Spistillat")

#RIGHT NOW THIS JUST DOES MILLEPORA COMPARISONS, BUT COULD ADD

for (compare.spp in spp.list){
  dnds = read.table(paste('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_', compare.spp, '.txt', sep = ""), header = T)#for digitifera
  dnds.dat = merge.dats(mdat, dnds, orthos, iso2seq, 'dN', F)
  head(dnds.dat)
  x = c()
  for (i in dnds.dat$EST){
    num = strsplit(i, "=")[[1]][2]
#     print(num)
    res = paste("isogroup", num, sep = "")
    x = append(x, res)
  }
  dnds.out = data.frame(x, dnds.dat$dnds)
  colnames(dnds.out) = c('EST', 'dNdS')
  head(dnds.out)
  write.table(dnds.out, paste("../kog/", compare.spp, "-Amillepora_dnds.txt", sep = ""), row.names = F, sep = ",", quote = F)
}

###################################################################
################# CORRELATION WITH dN WITHIN GO CATEGORIES ########
###################################################################

x = read.table('amil_iso2go.tab', header = F)
#THE ORIGINAL TABLE HAS THE GO LISTINGS, WHICH ARE CUMBERSOME, BUILD NEW FRAME WITH JUST ISOGROUPS AND PROCESS NAME
pdat = data.frame(x$V1, x$V3)
colnames(pdat) = c("EST", "process")
head(pdat)
bpdat = merge(pdat, dnds.dat, by = 'EST')
#GET THE PROCESS NAMES
ps = levels(pdat$process)
#SET UP EMPTY VECTORS TO STORE MEAN AND STANDARD ERROR IN FOR EACH BIOLOGICAL PROCESS
for (i in ps){
  print(i)
  sub = bpdat[bpdat$process == i,]
  if (nrow(sub) < 50){
    print("Not enough genes")
    next
  }
#   plot(dN~mbd.score, data = sub, main = i)
#   lm1 = lm(dN~mbd.score, data = sub)
#   print(summary(lm1))
#   abline(lm1, col = 'red')
  met.sub = sub[sub$mbd.score > 0,]
  u.sub = sub[sub$mbd.score < 0,]
  print(t.test(met.sub$dN, u.sub$dN, alternative = 'less'))
  mns = c(mean(met.sub$dN), mean(u.sub$dN))
  sterrs = c(std.error(met.sub$dN), std.error(u.sub$dN))
  xs = c(1,2)
  plotCI(xs, mns, uiw = sterrs, add = F, scol = 'grey', lwd = 1, main = i, xlim = c(0, 3), axes = F, ylab = 'dN')
  axis(1, at = c(1,2), labels = c('Met', 'UB'))
  axis(2)
}
##not very impressive but more or less the right trend

################################################################################
################ LOOKING AT CPGoe AND MUTATION RATES ###########################
################################################################################

#IMPORT DATA
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/")
#read in the species list
species = read.table("speciesList.txt")
colnames(species) = c('sp')
dat = data.frame()
#Import CpG data output from 
for (i in 1:length(species[,1])){
  sp = species[i,'sp']
  print(paste("Adding Data for", sp, "..."))
  #for each species read in the CpG data and append to dataframe
  x = read.table(paste(sp, "CpG.txt", sep = ""), header = T)
  x$sp = substr(sp, 1, 10)
  x$rep.strat = species[i,'rep.strat']
  dat = rbind(dat,x)
}
head(dat)

##############################################################
################ FILTER CPGoe DATA ###########################
##############################################################

#SET FILTERING PARAMETERS
minlength = 300
maxlength = 20000 ###original value for gene bodies was 20000
minOE = .001  ##.001
maxOE = 2
#filter the dataset
filter = function(cg, minlength, maxlength, maxOE, minOE, title){
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
  cg = cg[cg$cpgOE>minOE,] ##what are these?
  cg = na.omit(cg)
  return(cg)
}
cgdat= filter(dat, minlength, maxlength, maxOE, minOE, "All Species")

#PLOT DISTRIBUTION FOR ALL DATA
hist(cgdat$cpgOE, breaks = 50, main = "All Species", prob = T, ylim = c(0, 1.4), xlab = "CpGoe", xlim = c(0,1.5))

##PLOT THE DISTRIBUTION FOR EACH INDIVIDUAL SPECIES
for (i in 1:length(species[,1])){
  sp = substring(as.character(species[i,1]), 1, 10)
  print(sp)
  st = species[i,2]
  sub = cgdat[cgdat$sp == sp,]
  hist(sub$cpgOE, breaks = 50, prob = T, main = paste(sp, st, nrow(sub), 'genes'))
}
#OUTPUT THE DATASET
# write.table(cgm,'/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores7-27-14__1000LengthFilter.out',row.names=FALSE, quote = F)

##############################################################
################ PLOT RELATIONSHIP WITH CPGOE ################
##############################################################

merge.cpg.dats = function(compare.species, SPECIES){
  sub.compare = substring(compare.species, 1, 10)
  comparedat = cgdat[cgdat$sp == sub.compare,]
  head(comparedat)
  SUB.SPECIES = substring(SPECIES, 1, 10)
  #pull the pairwise substitution rates for a given species
  dnds = read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", SUB.SPECIES, '.txt', sep = ""), header = T)
  head(dnds)
  #subset the substition rates for just our compare species
  sub.dnds = dnds[dnds$species == sub.compare,]
  head(sub.dnds)
  head(orthos)
  sub.orth = merge(sub.dnds, orthos, by = 'EST')
  head(sub.orth)
  sub.orth2 = data.frame(sub.orth$EST, sub.orth$dN, sub.orth$dS, sub.orth[,sub.compare])
  colnames(sub.orth2) = c('dig', 'dN', 'dS', 'EST')
  head(sub.orth2)
  result = merge(sub.orth2, comparedat, by = 'EST')
  head(result)
  return(result)
}

plot.cpg.mut = function(compare.species, speciesList){
  compare.species <- substring(compare.species, 1, 10)
  for (sp in speciesList){
    sp <- substring(sp, 1, 10)
    if (sp == compare.species){
      next
    }
    dat = merge.cpg.dats(compare.species, sp)
    title <- paste(compare.species, 'vs', sp)
    plot.mut(dat, MUT.TYPE, 'cpgOE', MUT.TYPE, 'cpgOE', title, FALSE, 0.65)
  }
}

speciesList <- c('Amillepora', 'Adigitifera', 'Apalmata', 'Ahyacinthus', 'Atenuis', 'Ssiderea', 'Pastreoides', 'Mcavernosa', 'Pcarnosus', 'Pstrigosa', 'Fscutaria', 'Mauretenra', 'Pdamicornis', 'Spistillata', 'Shystrix')
speciesList2 <- c('Amillepora', 'Apalmata', 'Ahyacinthus', 'Atenuis', 'Ssiderea', 'Pastreoides', 'Mcavernosa', 'Pcarnosus', 'Pstrigosa', 'Fscutaria', 'Mauretenra', 'Pdamicornis', 'Spistillata', 'Shystrix')

quartz()
par(mfrow = c(3,5))
plot.cpg.mut('Apalmata', speciesList)
for (ANCHOR in speciesList2){
  quartz()
  par(mfrow = c(3,5))
  plot.cpg.mut(ANCHOR, speciesList)
}

merge.cpg.dats()

################ PLOT dN or dS FOR ANY COMPARISON SPECIES AND A CHOSEN FAMILY ######################
MUT.TYPE = 'dN'
FAMILY <- acroporids
FAMILY <- porites
FAMILY <- favids
COMPARE.SP = 'Ssiderea'

quartz()
for (sp in species$sp){
  if (sp == COMPARE.SP){
    next
  }
  FILENAME = paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/figures/Ssiderea/dS/", paste(sp, '.pdf', sep = ""), sep = "")
#   pdf(FILENAME)
  z<-plot.cpg.mut(COMPARE.SP, sp, TRUE)
#   dev.off()
}
 
plot.CPG.mut.stack = function(ddat, Ycol, Xcol, XLAB, MAIN, rev, SEPERATOR, LINE.COLOR, ITERATION, XLIM, YLIM){
  window_data = window_plot(1/8, Xcol, Ycol, ddat, Cutoff)
  if (ITERATION < 2){
    plot(mn ~ x, type = "n", xlim = XLIM, ylim = YLIM, data = window_data, axes = F, xlab = "CpGoe", ylab = YLAB)
    axis(1)
    axis(2)
    loess_fit <- loess(mn ~ x, window_data, span = .9, se = T)
    lines(window_data$x, predict(loess_fit), lwd=2, col = LINE.COLOR)
  }
  if (ITERATION >= 2){
    loess_fit <- loess(mn ~ x, window_data, span = .9, se = T)
    lines(window_data$x, predict(loess_fit), lwd=2.5, col = LINE.COLOR)
  }
}#function for reformatting and plotting data


z = read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", "Ssiderea", '.txt', sep = ""), header = T)
plot.CPG.mut.stack(z, MUT.TYPE, cpgOE, 'cpgoe', 'title', TRUE, .65, )





















###### PLOT THE CGM NONACROPORIDS ##################
speciesList <- c('Ssiderea', 'Pastreoide', 'Mcavernosa', 'Pcarnosus', 'Pstrigosa', 'Fscutaria', 'Mauretenra', 'Pdamicorni', 'Spistillat', 'Shystrix')
compare.species = 'Ssiderea'
sub.compare = substring(compare.species, 1, 10)
comparedat = cgdat[cgdat$sp == sub.compare,]
for (DUDE in speciesList){
  print(DUDE)
  sub.dude = substring(DUDE, 1, 10)
  #pull the pairwise substitution rates for this species
  dnds.dat = read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", sub.dude, '.txt', sep = ""), header = T)
  print(head(dnds.dat))
  #merge the substitution rate data with the cpg data for the compare species
  dnds.dat = merge.cgm(comparedat, dnds, orthos, iso2seq, 'dN', TRUE, DUDE)
  window.dat = plot.mut(dnds.dat, MUT.TYPE, 'cpgOE', MUT.TYPE, 'mbd.score', species, FALSE, 0.65)
}

#########Plotting Parameters 
quartz()
par(mfrow=c(1,1))
par(mar=c(5, 4, 4, 8) + 0.1)
par(xpd=F)
XLAB = "mbd.score"
YLAB = "dN"
FILT = T

###dN specific plotting parameters
MUT.TYPE = 'dN'
mut.ylim = c(0.004, 0.016)
YLAB = "dN"

###dS specific plotting parameters
MUT.TYPE = 'dS'
mut.ylim = c(.025, .125)
YLAB = "dS"

###PLOT THE LOESS LINES FOR THE SLIDING WINDOW FOR THE ACROPORIDS
speciesList <- c("Adigitifera", "Ahyacinthus", "Apalmata", "Atenuis")
par(xpd=F)
COLORS = rainbow(length(speciesList))
i = 0
for (species in speciesList){
  i = i + 1
  sub.species = substr(species, 1, 10)
  line.color = COLORS[i]
  print(species)
  dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", sub.species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
  plot.mut.stack(dnds.dat, MUT.TYPE, 'cpgOE', XLAB, species, TRUE, 0, line.color, i, c(2, -2), mut.ylim)
}
par(xpd=T)
legendY = mut.ylim[1] + (mut.ylim[2] - mut.ylim[1])*4/5
legend(-2.25, legendY, rev(speciesList), fill = rev(COLORS))






