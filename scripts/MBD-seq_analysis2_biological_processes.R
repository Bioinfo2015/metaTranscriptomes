#MBD-seq_analysis.R
#Groves Dixon
#1/16/15

#upload R object
load('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/MBD-seq_Image.R')
setwd('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/')
############################################################################
################## PLOT BIOLOGICAL PROCESSES ###############################
############################################################################
#Here we generate a plot for the mean MBD score for a series of biological processes
par(mfrow = c(1,1))
require(plotrix)
#READ IN THE DATA
x = read.table('process2iso.txt', header = F)# see appendix 'setting up go term data' in Methylation_Pulldown_Analysis_Walkthrough.txt
#THE ORIGINAL TABLE HAS THE GO LISTINGS, WHICH ARE CUMBERSOME, BUILD NEW FRAME WITH JUST ISOGROUPS AND PROCESS NAME
pdat = data.frame(x$V1, x$V3)
colnames(pdat) = c("isogroup", "process")
head(pdat)
#GET THE PROCESS NAMES
ps = levels(pdat$process)
#SET UP EMPTY VECTORS TO STORE MEAN AND STANDARD ERROR IN FOR EACH BIOLOGICAL PROCESS
means = c()
st.errs = c()
Ns = c()
for (i in ps){
  sub = pdat[pdat$process == i,]
  comb = merge(sub, mdat, by = 'isogroup')
  comb$mbd.score = comb$mbd.score*-1#this is because reversing the X axis screws up plotCI()
  means = append(means, mean(comb$mbd.score))
  st.errs = append(st.errs, std.error(comb$mbd.score))
  Ns = append(Ns, nrow(comb))
}#populates the vectors with the means, Ns and standard errors
y = data.frame(ps, means, st.errs, Ns)
y = y[with(y, order(means)),]
quartz()
par(mar=c(3,14,1,1))
plotCI(y$means, 1:length(y$means), uiw = y$st.errs, liw = y$st.errs, err="x", yaxt="n", pch=19, lwd=2, sfrac = .0075, axes = F, ylab = "")
axis(2, at = 1:nrow(y), labels = F, tick = T, xlim = c(.5, 1.3))
axis(1)
#Look at the GO term order
#Then type out nicely formated form
go.labels = c('DNA metabolism', 'Ribosome biogenesis', 'Translation', 'RNA metabolism', 'Transcription', 'Cell cycle and proliferation', "Stress response", "Transport", "Development", "Defense response", "Cell adhesion", "Signal transduction", "Response to stimulus", "Cell-cell signaling")
y.labels = paste(go.labels," (",y$Ns,")", sep = "")
mtext(y.labels, side = 2, line = .75, at = 1:nrow(y), outer=F,las=1, cex = 1.1)

#IN ADDITION TO PLOTTING THESE, WE CAN RUN MISHA'S GO ENRICHMENT TEST DIRECTLY ON THE LOG2FOLD CHANGES
#THIS REQUIRES THAT WE OUTPUT THE DATA IN A PARTICULAR FORMAT. THEN ANALYZE USING GO_MWU_MBD-seq_enrichment.R
out = data.frame(mdat$EST, mdat$log2FoldChange)
colnames(out) = c('gene', 'log2FoldChange')
head(out)
corrected.iso = c()
for (i in 1:length(out$gene)){
  
  corrected.iso = append(corrected.iso, paste(strsplit(as.character(out$gene[i]), "=")[[1]][1], strsplit(as.character(out$gene[i]), "=")[[1]][2], sep = ""))
}
out$gene = corrected.iso
write.table(out, "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/go/MBD_foldChangeForGO_V2_5-14-15.txt", sep = ',', quote = F, row.names = F)
#write one for KOGS to
write.table(out, "/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/kog/MBD_foldChangeForKOG_V2_5-14-15.txt", sep = ',', quote = F, row.names = F)
