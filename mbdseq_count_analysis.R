#ANALYZING MBD-SEQ MAPPING COUNTS
#THIS ASSUMES YOU HAVE MAPPED THE READS AGAINST A TRANSCRIPTOME
#AND RUN getCpGoe.py ON THE TRANSCRIPTOME

#UPLOAD THE COUNTS DATA
setwd("/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/pull_down")
dat = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/pull_down/allcounts.txt', header = T)
colnames(dat) = c("met1", "ub1", "met2", "ub2")#note this is correcting for the mislabeling before sequencing

#RANK EACH GENE FOR EACH REPLICATE BASED ON COUNTS
dat$met1.r = rank(dat$met1)
dat$met2.r = rank(dat$met2)
dat$ub1.r = rank(dat$ub1)
dat$ub2.r = rank(dat$ub2)
dat$EST = rownames(dat)

#NORMALIZE THE COUNTS FOR EACH REPLICATE BY DIVIDING THE COUNTS BY THE REPLICATE'S MEAN COUNTS
normalize = function(vector){
  norm = vector/mean(vector)
  return(norm)
}

#GET THE AVERAGE OF THE NORMALIZED COUNTS BETWEEN THE TWO REPLICATES FOR EACH TREATMENT
#THEN TURN THE MEAN VALUES INTO RANKS
dat.n = data.frame(normalize(dat$met1), normalize(dat$met2), normalize(dat$ub1), normalize(dat$ub2))
dat.n$mnMet = apply(dat.n[,1:2], 1, mean)
dat.n$mnUb = apply(dat.n[,3:4], 1, mean)
dat.n$Rmet = rank(dat.n$mnMet)
dat.n$Rub = rank(dat.n$mnUb)
dat.n$EST = dat$EST
colnames(dat.n) = c('met1','met2','ub1','ub2','mnMet','mnUb', 'Rmet','Rub','EST')
head(dat.n)

#UPLOAD THE CPG DATA. THIS WAS OUTPUT FROM CpG_distribution.R
cgm = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores8-22-14__CDS.txt", header = T)
head(cgm)

#MERGE THE DATAFRAMES SO THAT YOU HAVE CpGoe DATA ALONGSIDE THE COUNTS DATA
dat2 = merge(dat.n, cgm, by = "EST")

#GET THE DIFFERENCE BETWEEN THE PULLDOWN RANKS AND THE UNBOUND RANKS
#THIS IS OUR MEANSURE OF METHYLATION STATE
dat2$mbd.score = dat2$Rmet - dat2$Rub

#SET UP FUNCTION TO PLOT LINEAR RELATIONSHIPS
linear = function(y, x, dat){
  plot(dat[,x]~dat[,y], main = x)
  lm1 = lm(dat[,x]~dat[,y])
  abline(lm1, col = 'red')
  print(summary(lm1))
}

#PLOT CORRELATION WITH CPGOE
quartz()
par(mfrow = c(1,2))
plot(mbd.score~cpgOE, data = dat2, pch = 19, cex = 0.1, ylab = "Mean Read Count Rank mbd.scoreerence (pulldown - flowthru)", main = "Correlation with CpGoe")
plot(mbd.score~gpcOE, data = dat2, pch = 19, cex = 0.1, ylab = "Mean Read Count Rank mbd.scoreerence (pulldown - flowthru)", main = "Correlation with GpCoe (control)")
lm1 = lm(mbd.score~cpgOE, data = dat2)
summary(lm1)
abline(lm1, col = 'red', lwd = 2)

#PLOT THE DENSITY PLOT FOR THE CORRELATION
require(ggplot2)
qplot(mbd.score, cpgOE, data = dat2, geom = "density2d")
g = ggplot(dat2, aes(mbd.score, cpgOE))
g + geom_point() + geom_density2d()
g + stat_density2d(geom="tile", aes(fill = ..density..), contour = F) + scale_fill_continuous(low = 'blue', high = 'green', guide = 'colorbar')# + coord_cartesian(xlim = c(0, 1.35), ylim = c(-20000,20000))

#REVERSE THE AXES SO IT LOOKS PRETTIER
g = ggplot(dat2, aes(cpgOE, mbd.score))
g + geom_point(xlab = 'hey') + geom_density2d()
g + stat_density2d(geom="tile", aes(fill = ..density..), contour = F) + scale_fill_continuous(low = 'blue', high = 'green', guide = 'colorbar') + coord_cartesian(xlim = c(0, 1.35), ylim = c(-20000,20000))

#OUTPUT THE DATAFRAME
write.table(dat2, "/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/pull_down/mbdScoreIsogroups.txt", quote = F, row.names = F)
