#MBD-seq_analysis4_dNdS.R


#upload the data
setwd('/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory')
load("MBD-seq_Image.R")
source("~/git_Repositories/metaTranscriptomes/scripts/MBD-seq_analysis_functions.R")
dnds = read.table('pair-wise_dNdS_Amillepora.txt', header = T)#for digitifera
colnames(dnds) = c('contig', 'species', 'dN', 'dS', 'dNdS')

#Get isogroup names for the contigs
head(iso2seq)
head(dnds)
sub = merge(dnds, iso2seq, by = 'contig')
head(sub)
dnds.dat = merge(sub, mdat, by = 'isogroup')[,1:7]
head(dnds.dat)

#make a preliminary plot for Adigitifera
require(plotrix)
QUANTILE = TRUE
Cutoff = 50
mut.type = 'dN'
meth.type = 'mbd.score'
dig = filter.sub.rates(dnds.dat, 'Adigitifer', mut.type, meth.type, TRUE, TRUE, 2, 'dS')
wdat = window_plot(1/15, 'mbd.score', 'dN', dig, 10)
plotCI(wdat$x, wdat$mn, uiw = wdat$sterr, liw = wdat$sterr, add = F, scol = 'black', lwd = 1, main = paste('dig', nrow(sub)), xlim = c(2,-2), pch = 26)

#choose global plotting variables to plot
#pair-wise substitution rates vs methylation scores
mut.type = 'dN'
filter.mut = 'dS'
meth.type = 'mbd.score'    #mbd.score or cpg
filter = T
remove.zeros = T
factor = 2
span = .5
loess = T
window = 1/20
line.width = 1.5
cex.axis = 1.25
cex.lab = 1.25
title.line = .6

#plot for acroporids
quartz()
par(mfrow = c(2,2))
par(mar = c(4,3.5,1.5,.5) + 0.1)
acro = c("Adigitifer", "Ahyacinthu", "Apalmata", "Atenuis")
acro.right = c("A.digitifera", "A.hyacinthus", "A.palmata", "A.tenuis")
# speciesList = acro
# for (i in 1:length(speciesList)){
  # sp = as.character(speciesList[i])
  # sub = filter.sub.rates(dnds.dat, sp, mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
  # sp = as.character(acro.right)[i]
  # win = do.window.plot(window, meth.type, mut.type, sub, 10, "\n", "\n", '\n', loess, span, limits = F)
  # title(paste(sp, nrow(sub)), font.main = 3, line = .75)
  # title(xlab = "MBD-score", ylab = mut.type, cex.lab = cex.lab, line = 2.25)
  # yat = signif(seq(from = min(win$mn - win$sterr), to = max(win$mn + win$sterr), by = (max(win$mn + win$sterr) - min(win$mn - win$sterr))/2), 2)
  # axis(1, cex.axis = cex.axis, at = c(-2, 0, 2))
  # axis(2, cex.axis = cex.axis, at = yat)
  # # axis(2)
# }

#set margins for top figures
quartz()
par(mfrow = c(2,2))
par(mar = c(2.5,3.5,1.5,.6) + 0.1)
#Adig
line.width = 2
sub = filter.sub.rates(dnds.dat,'Adigitifer', mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
sp = "A.digitifera"
win = do.window.plot(window, meth.type, mut.type, sub, 10, "\n", "\n", '\n', loess, span, limits = F)
title(paste(sp, nrow(sub)), font.main = 3, line = title.line)
title(xlab = "\n", ylab = mut.type, cex.lab = cex.lab, line = 2.25)
yat = signif(seq(from = min(win$mn - win$sterr), to = max(win$mn + win$sterr), by = (max(win$mn + win$sterr) - min(win$mn - win$sterr))/2), 2)
axis(1, cex.axis = cex.axis, at = c(-2, 0, 2))
axis(2, cex.axis = cex.axis, at = yat)
line.width = 3
plot.meth.bars(sub, 'mbd.score', 'dN', separator, c('green', 'red'))
#Ahyacinthus
line.width = 2
sub = filter.sub.rates(dnds.dat,'Ahyacinthu', mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
sp = "A.hyacinthus"
win = do.window.plot(window, meth.type, mut.type, sub, 10, "\n", "\n", '\n', loess, span, limits = F)
title(paste(sp, nrow(sub)), font.main = 3, line = title.line)
title(xlab = "\n", ylab ="\n", cex.lab = cex.lab, line = 2.25)
yat = signif(seq(from = min(win$mn - win$sterr), to = max(win$mn + win$sterr), by = (max(win$mn + win$sterr) - min(win$mn - win$sterr))/2), 2)
axis(1, cex.axis = cex.axis, at = c(-2, 0, 2))
axis(2, cex.axis = cex.axis, at = yat)
line.width = 3
plot.meth.bars(sub, 'mbd.score', 'dN', separator, c('green', 'red'))

par(mar = c(4,3.5,1.5,.5) + 0.1) #set margins for bottom figures
#Apalmata
line.width = 2
sub = filter.sub.rates(dnds.dat,'Apalmata', mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
sp = "A.palmata"
win = do.window.plot(window, meth.type, mut.type, sub, 10, "\n", "\n", '\n', loess, span, limits = F)
title(paste(sp, nrow(sub)), font.main = 3, line = title.line)
title(xlab = "MBD-score", ylab ="dN", cex.lab = cex.lab, line = 2.25)
yat = signif(seq(from = min(win$mn - win$sterr), to = max(win$mn + win$sterr), by = (max(win$mn + win$sterr) - min(win$mn - win$sterr))/2), 2)
axis(1, cex.axis = cex.axis, at = c(-2, 0, 2))
axis(2, cex.axis = cex.axis, at = yat)
line.width = 3
plot.meth.bars(sub, 'mbd.score', 'dN', separator, c('green', 'red'))
#Atenuis
line.width = 2
sub = filter.sub.rates(dnds.dat,'Atenuis', mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
sp = "A.tenuis"
win = do.window.plot(window, meth.type, mut.type, sub, 10, "\n", "\n", '\n', loess, span, limits = F)
title(paste(sp, nrow(sub)), font.main = 3, line = title.line)
title(xlab = "MBD-score", ylab ="\n", cex.lab = cex.lab, line = 2.25)
yat = signif(seq(from = min(win$mn - win$sterr), to = max(win$mn + win$sterr), by = (max(win$mn + win$sterr) - min(win$mn - win$sterr))/2), 2)
axis(1, cex.axis = cex.axis, at = c(-2, 0, 2))
axis(2, cex.axis = cex.axis, at = yat)
line.width = 3
plot.meth.bars(sub, 'mbd.score', 'dN', separator, c('green', 'red'))



#plot for non acroporid corals
mut.type = "dN"
window = 1/20
line.width = 2
non.acro = c('Ssiderea', 'Pastreoide', 'Fscutaria', 'Mcavernosa', 'Pstrigosa', 'Mfaveolata', 'Pdaedalea', 'Pcarnosus', 'Mauretenra', 'Pdamicorni', 'Spistillat', 'Shystrix')
format.list = c('S.siderea', 'P.astreoides', 'F.scutaria', 'M.cavernosa', 'P.strigosa', 'O.faveolata', 'P.daedalea', 'P.carnosus', 'M.auretenra', 'P.damicornis', 'S.pistillata', 'S.hystrix')
quartz()
par(mfrow = c(3,4))
speciesList = non.acro
percent.dif.list = c()
for (i in 1:length(speciesList)){
  sp = as.character(speciesList[i])
  sub = filter.sub.rates(dnds.dat, sp, mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
  print(head(sub))
  sp = format.list[i]
  do.window.plot(window, meth.type, mut.type, sub, 10, meth.type, mut.type, paste(sp, nrow(sub)), T, span, limits = F)
  axis(1)
  axis(2)
  line.width = 3
  pct = plot.meth.bars(sub, 'mbd.score', mut.type, separator, c('green', 'red'))
  line.width = 2
  percent.dif.list = append(percent.dif.list, pct)
}
#look at the mean percent increase between categories
mean(percent.dif.list); std.error(percent.dif.list)

#plot anemones
quartz()
par(mfrow = c(1,3))
actin = c('Nvectensis', 'Apallida', 'Aelegantis')
format.list = c('N.vectensis', 'A.pallida', 'A.elegantissima')
percent.dif.list = c()
speciesList = actin
line.width = 2
for (i in 1:length(speciesList)){
  sp = as.character(speciesList[i])
  sub = filter.sub.rates(dnds.dat, sp, mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
  print(head(sub))
  sp = format.list[i]
  do.window.plot(window, meth.type, mut.type, sub, 10, meth.type, mut.type, paste(sp, nrow(sub)), T, span, limits = F)
  axis(1)
  axis(2)
  line.width = 3
  pct = plot.meth.bars(sub, 'mbd.score', mut.type, separator, c('green', 'red'))
  line.width = 2
  percent.dif.list = append(percent.dif.list, pct)
}
#look at the mean percent increase between categories
mean(percent.dif.list); std.error(percent.dif.list)

#build plots with lines stacked
#dN
mut.type = 'dN'
filter.mut = 'dS'
window = 1/10
line.width = 3
remove.zeros = T
filter = T
cex.X = 2
cex.Y = 1.75
xlim = c(-2.6, 3)
colors = rainbow(24)
par(mfrow = c(3, 1))
par(mar = c(2, 4, 1, -.1) + 0.1)
plot.stacked.lines(acro, colors[20:24], xlim, c(.0048, .0265))
plot.stacked.lines(non.acro, colors[6:18], xlim, c(.10, .235))
plot.stacked.lines(actin, colors[1:3], xlim, c(.25, .41))

#dS
quartz()
par(mfrow=c(3,1))
mut.type = 'dS'
plot.stacked.lines(acro, colors[20:24], xlim, c(.025, .1))
plot.stacked.lines(non.acro, colors[6:18], xlim, c(1.2, 3.0))
plot.stacked.lines(actin, colors[1:3], xlim, c(23, 37))

#dNdS
quartz()
par(mfrow=c(3,1))
mut.type = 'dNdS'
filter.mut = 'dNdS'
plot.stacked.lines(acro, colors[20:24], xlim, c(.17, .326))
plot.stacked.lines(non.acro, colors[6:18], xlim,  c(0.058, 0.1))
plot.stacked.lines(actin, colors[1:3], xlim,  c(0.021, 0.033))


#plot the pallets for color matching later
plot.pallet(acro, colors[20:24])
plot.pallet(non.acro, colors[6:18])
plot.pallet(actin, colors[1:3])




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
      plot(mn ~ x, type = "n", xlim = XLIM, ylim = mut.ylim, data = window_data, axes = F, xlab = "log2FoldChange", ylab = YLAB)
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
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dS'
MUT.YLIM = c(.025, .125)
YLAB = "dS"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dNdS
MUT.TYPE = 'dnds'
MUT.YLIM = NULL
YLAB = "dN/dS"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

############### PLOT FOR OTHER CORALS ###############
#dN
###PLOT dN FOR NON ACROPORID SPECIES
quartz()
par(mfrow = c(1, 3))
COLS = COLOR.SET[5:8]
MUT.TYPE = 'dN'
MUT.YLIM = NULL
YLAB = "dN"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dS'
MUT.YLIM = c(1.8, 3.2)
YLAB = "dS"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dN/dS
MUT.TYPE = 'dnds'
MUT.YLIM = c(0.048, 0.09)
YLAB = "dN/dS"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

########## PLOT dS FOR ACTINIARIA #################
#dN
COLS = COLOR.SET[9:10]
quartz()
par(mar=c(5, 4, 0, 8) + 0.1)
par(xpd=F)
MUT.TYPE = 'dN'
MUT.YLIM = c(.185, .235)
YLAB = "dN"
plot.mut.stackF(actin, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
quartz()
par(mar=c(5, 4, 0, 8) + 0.1)
par(xpd=F)
MUT.TYPE = 'dS'
MUT.YLIM = c(24.5, 30)
YLAB = "dS"
plot.mut.stackF(actin, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

#PLOT THE LOESS LINES FOR PAIR-WISE MUTATION RATES FOR THE ACROPORID SPECIES ########

#########Plotting Parameters 
quartz()
par(mfrow=c(1,1))
par(mar=c(5, 4, 4, 8) + 0.1)
par(xpd=F)
XLAB = "log2FoldChange"
acroporids <- c("Adigitifera", "Ahyacinthus", "Apalmata", "Atenuis")
FILT = T

###dN specific plotting parameters
MUT.TYPE = 'dN'
MUT.YLIM = c(0.004, 0.016)
YLAB = "dN"
plot.mut.stackF(acroporids, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2))

###dS specific plotting parameters
MUT.TYPE = 'dS'
MUT.YLIM = c(.025, .125)
YLAB = "dS"
plot.mut.stackF(acroporids, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2))

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
  #   plot(dN~log2FoldChange, data = sub, main = i)
  #   lm1 = lm(dN~log2FoldChange, data = sub)
  #   print(summary(lm1))
  #   abline(lm1, col = 'red')
  met.sub = sub[sub$log2FoldChange > 0,]
  u.sub = sub[sub$log2FoldChange < 0,]
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

