#MBD-seq_analysis2_biological_processes.R
#This script plots the mean and standard error MBD-score for a set of 
#biological process GO terms 
#Groves Dixon
#1/16/15

#upload R object
setwd('/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory')
load('MBD-seq_Image.R')
require(plotrix)

#read in the data
x = read.table('process2iso.txt', header = F) #this is a table linking each isogroup with a set of biological process GO terms

#grab subset of original bulky dataframe
pdat = data.frame(x$V1, x$V3)
colnames(pdat) = c("isogroup", "process")
head(pdat)

#assign the process names
ps = levels(pdat$process)

#gather the mean, standard error, and number of genes for each process
means = c()
st.errs = c()
Ns = c()
for (i in ps){
  sub = pdat[pdat$process == i,]
  comb = merge(sub, mdat, by = 'isogroup')
  comb$mbd.score = comb$mbd.score
  means = append(means, mean(comb$mbd.score))
  st.errs = append(st.errs, std.error(comb$mbd.score))
  Ns = append(Ns, nrow(comb))
}#populates the vectors with the means, Ns and standard errors

#assemble into a dataframe
y = data.frame(ps, means, st.errs, Ns)
y = y[with(y, order(means)),]

#plot the results
quartz()
par(mar=c(5,14,1,1))
plotCI(y$means, 1:length(y$means), uiw = y$st.errs, liw = y$st.errs, err="x", yaxt="n", pch=26, lwd=3, sfrac = .0075, axes = F, ylab = "", xlab = "\n")
axis(2, at = 1:nrow(y), labels = F, tick = T, xlim = c(.5, 1.3))
axis(1, cex.axis = 1.25)
title(xlab = "MBD-score", cex.lab = 1.25)
#Look at the GO term order
#Then type out nicely formated form
go.labels = rev(c('DNA metabolism', 'Ribosome biogenesis', 'Translation', 'RNA metabolism', 'Transcription', 'Cell cycle and proliferation', "Stress response", "Transport", "Development", "Defense response", "Cell adhesion", "Signal transduction", "Response to stimulus", "Cell-cell signaling"))
y.labels = paste(go.labels," (",y$Ns,")", sep = "")
mtext(y.labels, side = 2, line = .75, at = 1:nrow(y), outer=F,las=1, cex = 1.1)

