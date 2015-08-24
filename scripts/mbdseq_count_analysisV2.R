#ANALYZING MBD-SEQ MAPPING COUNTS
#THIS ASSUMES YOU HAVE MAPPED THE READS AGAINST A TRANSCRIPTOME
#AND RUN getCpGoe.py ON THE TRANSCRIPTOME

#UPLOAD THE COUNTS DATA
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/")
dat = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/MBD-seq/data_files/counts_2-10-15_dupsIncluded.txt', header = T)
colnames(dat) = c("met1", "met2", "ub1", "ub2")#note this is correcting for the mislabeling before sequencing
head(dat)

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
dat.n$diff = dat.n$mnMet - dat.n$mnUb
dat.n$Rmet = rank(dat.n$mnMet)
dat.n$Rub = rank(dat.n$mnUb)
dat.n$Rdiff = rank(dat.n$diff)
dat.n$EST = dat$EST
colnames(dat.n) = c('met1','met2','ub1','ub2','mnMet','mnUb','mnDiffs', 'Rmet','Rub','RmnDiffs', 'EST')
head(dat.n)

#UPLOAD THE CPG DATA. THIS WAS OUTPUT FROM CpG_distribution.R
cgm = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores8-22-14__CDS.txt", header = T)
# cgm = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores1-13-15__CDSnosub.txt", header = T)##########this one was to look at the nonsubsetted one, but the subset looks way better

#MERGE THE DATAFRAMES SO THAT YOU HAVE CpGoe DATA ALONGSIDE THE COUNTS DATA
dat2 = merge(dat.n, cgm, by = "EST")

#GET THE DIFFERENCE BETWEEN THE PULLDOWN RANKS AND THE UNBOUND RANKS
#THIS IS OUR MEANSURE OF METHYLATION STATE
dat2$mbd.score = dat2$Rmet - dat2$Rub
dat2$mbd.score2 = dat2$mnMet / dat2$CpG
dat2$sanity.score2 = dat2$mnUb / dat2$CpG

##LOOK AT THE HISTOGRAMS
hist(dat2$mbd.score)
hist(log(dat2$mbd.score2))
hist(log(dat2$sanity.score2))

##LOOK AT INDIVIDUAL HISTOGRAMS
quartz()
par(mfrow = c(2,2))
hist(log(dat2$met1/dat2$CpG, 10), main = "Pulldown1 Distribution", xlab = 'log(Reads/CpG)', breaks = 20)
hist(log(dat2$met2/dat2$CpG, 10), main = "Pulldown2 Distribution", xlab = 'log(Reads/CpG)', breaks = 20)
hist(log(dat2$ub1/dat2$CpG, 10), main = "FlowThru1 Distribution", xlab = 'log(Reads/CpG)', breaks = 20)
hist(log(dat2$ub2/dat2$CpG, 10),  main = "FlowThru2 Distribution", xlab = 'log(Reads/CpG)', breaks = 20)




#SET UP FUNCTION TO PLOT LINEAR RELATIONSHIPS
linear = function(y, x, dat){
  plot(dat[,x]~dat[,y], main = x)
  lm1 = lm(dat[,x]~dat[,y])
  abline(lm1, col = 'red')
  print(summary(lm1))
}

#PLOT CORRELATION OF MET - UB RANKS CPGOE
quartz()
par(mfrow = c(1,2))
plot(mbd.score~cpgOE, data = dat2, pch = 19, cex = 0.1, ylab = "Mean Read Count Rank mbd.scoreerence (pulldown - flowthru)", main = "Correlation with CpGoe")
lm1 = lm(mbd.score~cpgOE, data = dat2)
summary(lm1)
abline(lm1, col = 'red', lwd = 2)
plot(mbd.score~gpcOE, data = dat2, pch = 19, cex = 0.1, ylab = "Mean Read Count Rank mbd.scoreerence (pulldown - flowthru)", main = "Correlation with GpCoe (control)")

#PLOT THE SECOND SCORING METHOD (mnMet / CpG counts)
quartz()
par(mfrow = c(1,2))
plot(log(mbd.score2, 10)~cpgOE, data = dat2, pch = 19, cex = 0.1, ylab = "log(Reads/CpG)", main = "Pulldown Correlation with CpGoe")
datLin = data.frame(dat2$mbd.score2, dat2$cpgOE)
colnames(datLin) = c('mbd.score2', 'cpgOE')
datLin = datLin[datLin$mbd.score2 > 0,]
lm1 = lm(log(datLin$mbd.score2, 10) ~ datLin$cpgOE)
summary(lm1)
# abline(lm1, col = 'red')
plot(log(sanity.score2, 10)~cpgOE, data = dat2, pch = 19, cex = 0.1, ylab = "log(Reads/CpG)", main = "Flowthru Correlation with CpGoe")
datLin = data.frame(dat2$sanity.score, dat2$cpgOE)
colnames(datLin) = c('mbd.score2', 'cpgOE')
datLin = datLin[datLin$mbd.score2 > 0,]
lm1 = lm(log(datLin$mbd.score2, 10) ~ datLin$cpgOE)
summary(lm1)
# abline(lm1, col = 'red')

#PLOT THE DENSITY PLOT FOR THE CORRELATION
require(ggplot2)
qplot(log(mbd.score2,10), cpgOE, data = dat2, geom = "density2d")
g = ggplot(dat2, aes(cpgOE, log(mbd.score2,10)))
g + geom_point() + geom_density2d()
g + stat_density2d(geom="tile", aes(fill = ..density..), contour = F) + scale_fill_continuous(low = 'black', high = 'orange', guide = 'colorbar')# + coord_cartesian(xlim = c(0, 1.35), ylim = c(-20000,20000))
#PLOT THE DENSITY PLOT SANITY
g = ggplot(dat2, aes(cpgOE, log(sanity.score2,10)))
g + geom_point() + geom_density2d()
g + stat_density2d(geom="tile", aes(fill = ..density..), contour = F) + scale_fill_continuous(low = 'black', high = 'orange', guide = 'colorbar')# + coord_cartesian(xlim = c(0, 1.35), ylim = c(-20000,20000))

#REVERSE THE AXES SO IT LOOKS PRETTIER
g = ggplot(dat2, aes(cpgOE, mbd.score))
g + geom_point(xlab = 'hey') + geom_density2d()
g + stat_density2d(geom="tile", aes(fill = ..density..), contour = F) + scale_fill_continuous(low = 'blue', high = 'green', guide = 'colorbar') + coord_cartesian(xlim = c(0, 1.35), ylim = c(-20000,20000))

#OUTPUT THE DATAFRAME
write.table(dat2, "/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/pull_down/mbdScoreIsogroups.txt", quote = F, row.names = F)
#===============================================================
head(dat2)



##########################################
#CREATE A LOG TRANSFORMED COLUMN FOR THE METHYLATION SCORE
#THIS IS TRICKY, CUZ SOME VALUES ARE ZERO, SO CHANGE THESE
#ASSIGN ALL ZEROS IN THE METHYLATION COUNTS TO THE MINIMUM NONZERO VALUE
#FUNCTION TO DO THIS:
log.trans = function(dat, col){
  x = dat[,col]/dat[,'cpgOE']
  z = x[x > 0]
  x = x + min(z)
  logtrans <- log(x, 10)
  x <- logtrans
  logtrans = logtrans - min(x)
  print(paste('This should be zero:', min(logtrans)))
  return(logtrans)
}
#EXECUTE THE FUNCTION FOR SCORE AND CONTROL SCORE
tdat$logscore = log.trans(tdat, 'mbd.score2')
tdat$logsanity= log.trans(tdat, 'sanity.score2')
tdat$logm1= log.trans(tdat, 'met1')
tdat$logm2= log.trans(tdat, 'met2')
tdat$logu1= log.trans(tdat, 'ub1')
tdat$logu2= log.trans(tdat, 'ub2')
######################################## 


##### TESTING FOR CONSISTENCIES WITH PREVIOUS DATA
head(dat2)
tdat = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/data_from_transplant_experiment/VSD_lbRT_nodup_nomiso__reduced.csv", header =T)
head(tdat)
tdat = merge(tdat, dat2, by = "EST")
length(tdat$EST)
##### SLIDING WINDOW FUNCTION #########
require(plotrix)
expression_plot = function(size, X, dat, cut){
  windows = quantile(dat[,X], probs = seq(0, 1, by = size))
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
    mn = append(mn,mean(sub$grandmeans))
    sterr = append(sterr,std.error(sub$grandmeans))
    x = append(x,windows[i])
    N = append(N, n)
  }
  print(paste(cut.count, "Windows were cut Because they had below", cut, "genes"))
  plot_dat = data.frame(mn,x,sterr, N)
  return(plot_dat)
}##function to get the window data to plot th
Cutoff = 10
window_data = expression_plot(.04, 'logu2', tdat, Cutoff) ## data for Carly's RNA seq
window_data = na.omit(window_data)##remove the windows that had no genes in them
plot(mn~x,data=window_data, main="", pch = 1, cex = 1, axes = T, cex.lab = 1)
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'grey', lwd = 0.5)

dat = tdat
X = 'logscore'
left = 0
right = 0.01603854
sub = dat[dat[,X]>=left,]
sub2 = sub[sub[,X]<right,]
n = length(sub2[,1])
n
head(tdat$logscore)
head(dat[,X])

nrow(tdat)


mn = round(mean(window_data$mn), digits = 1)
axis(1, at = c(0.2, 0.6, 1.0), cex.axis = CEX, mgp = c(1, .4 , 0))
mtext(expression("CpG"["O/E"]), side=1, line=1.6, outer=F, cex= CEX, font=1, adj = .48)
AT = c(4.3, 4.4, 4.5, 4.6)
axis(2, at = AT, label = T, tick = T, las = 1, cex.axis = CEX, mgp = c(1.8, .6, 0))
#mtext("Gene Expression", side=2, line=1.7, outer=F, cex= CEX, font=1, adj = .5)
loess_fit <- loess(mn ~ x, window_data, span = .8, se = T)
lines(window_data$x, predict(loess_fit),col="red",lwd=1)
green = mmix2$mean[1]
red = mmix2$mean[2]
int = separator
xs = c(green, int, red)
Y = YLIM[1] - 0.012
ys = c(Y, Y, Y)
colors = c("yellowgreen", "black", "red")
points(xs, ys, pch = 17, col = colors, cex = 2)









expression_plot = function(size, X, Y, dat, cut){
  windows = quantile(dat[,X], probs = seq(0, 1, by = size))
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
    x = append(x,windows[i])
    N = append(N, n)
  }
  print(paste(cut.count, "Windows were cut Because they had below", cut, "genes"))
  plot_dat = data.frame(mn,x,sterr, N)
  return(plot_dat)
}##function to get the window data to plot th
Cutoff = 10
window_data = expression_plot(.04, 'grandmeans', 'logu1', tdat, Cutoff) ## data for Carly's RNA seq
window_data = na.omit(window_data)##remove the windows that had no genes in them
plot(mn~x,data=window_data, main="", pch = 1, cex = 1, axes = T, cex.lab = 1)
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, scol = 'grey', lwd = 0.5)

head(dat2)
x = dat2$mbd.score2 - dat2$sanity.score2
hist(log(x, 10))
hist(tdat$logscore, breaks = 30)
hist(tdat$logsanity, breaks = 30)
