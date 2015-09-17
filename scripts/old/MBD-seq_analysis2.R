#################################################################################################################
################ DATA AND FUNCTIONS TO COMPARE PAIR-WISE SUBSTITUTION RATE VARIATION FOR A. MILLEPORA ###########
#################################################################################################################
#UPLOAD THE DNDS DATA FROM VARIOIUS OPTIONS
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/")
library(plotrix)
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
#   print(head(sub2))
  sub3 = data.frame(sub2$EST, sub2$Amillepora, sub2$dN, sub2$dS, sub2$dNdS)
  colnames(sub3) = c("digEST", "EST", "dN", "dS", "dnds")
  sub4 = merge(iso2seq, sub3, by = 'EST')
  colnames(sub4) = c('isotig', 'EST', 'digEST', 'dN', 'dS', 'dnds')
#   head(sub4)
  mut.dat = merge(met.dat, sub4, by = "EST")
  result = mut.dat
  print(paste("Number of Genes before Outlier Removal = ", nrow(result)))
  if (remove.outliers == TRUE){
    result = mut.dat[mut.dat[,mut.column] < (mean(mut.dat[,mut.column]) + SDs*sd(mut.dat[,mut.column])), ]
    outies = mut.dat[mut.dat[,mut.column] >=  (mean(mut.dat[,mut.column]) + SDs*sd(mut.dat[,mut.column])), ]
    assign('OUTLIERS', outies, envir = .GlobalEnv)
  }
  print(paste("Number of Genes after Outlier Removal = ", nrow(result)))
  return(result)
}
# SDs = .25
dnds.dat = merge.dats(mdat, dnds, orthos, iso2seq, 'dN', T)

#USE THE REFORMATTED DATA TO PLOT RELATIONSHIP BETWEEN PULLDOWN SCORE AND SUBSTITUTION RATES ##########
#FUNCTION TO PLOT A SLIDING WINDOW FOR GIVEN Y COLUMN AND X COLUMN
plot.mut = function(ddat, Ycol, Xcol, XLAB, YLAB, MAIN, rev, SEPERATOR){
#   ddat = ddat[abs(ddat[,Ycol]) != Inf,] #get rid of instances where dNdS was Inf because dS was 0
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
#     print(head(dnds.dat))
    N = nrow(dnds.dat)
    title = paste(species, "N =", N)
    window.dat = plot.mut(dnds.dat, MUT.TYPE, METH.TYPE, MUT.TYPE, METH.TYPE, title, TRUE, 0)
  }
}
plot_list_of_speciesCpGoe = function(speciesList){
  for (species in speciesList){
    print(species)
    dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
#     print(head(dnds.dat))
    dnds.dat = merge(cgm, dnds.dat, by = 'EST')
    N = nrow(dnds.dat)
    title = paste(species, "N =", N)
    window.dat = plot.mut(dnds.dat, MUT.TYPE, METH.TYPE, MUT.TYPE, METH.TYPE, title, FALSE, .65)
  }
}
#---------------------------------------------------------------------------------
#PLOT ACROPORIDS
METH.TYPE = 'log2FoldChange'
SDs = 2
MUT.TYPE = 'dN'
FILT = T
window.fraction = 1/10
quartz()
par(mfrow = c(1,4))
acro = c("Adigitifer", "Ahyacinthu", "Apalmata", "Atenuis")
METH.TYPE = 'log2FoldChange'
plot_list_of_species(acro)
METH.TYPE = 'cpgOE'
plot_list_of_speciesCpGoe(acro)

#PLOT OTHER CORALS
quartz()
par(mfrow = c(1,5))
non.acro = c('Ssiderea', 'Pastreoide', 'Pcarnosus', 'Spistillat', 'Pdamicorni')
METH.TYPE = 'log2FoldChange'
plot_list_of_species(non.acro)
METH.TYPE = 'cpgOE'
plot_list_of_speciesCpGoe(non.acro)
#PLOT ANOMONES
quartz()
par(mfrow = c(1,4))
actin = c('Nvectensis', 'Apallida')
METH.TYPE = 'log2FoldChange'
plot_list_of_species(actin)
METH.TYPE = 'cpgOE'
plot_list_of_speciesCpGoe(actin)

########## RUN FOR dS ##############

SDs = 2
MUT.TYPE = 'dS'
METH.TYPE = 'log2FoldChange'
quartz()
par(mfrow = c(1,4))
#PLOT ACROPORIDS
plot_list_of_species(acro)
#PLOT OTHER CORALS
quartz()
par(mfrow = c(1,5))
plot_list_of_species(non.acro)
#PLOT ANOMONES
quartz()
par(mfrow = c(1,4))
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
  i = 0
  percentages = c()
  for (species in speciesList){
    print(species)
    i = i + 1
    line.color = COLORS[i]
    dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
    dnds.dat2 = na.omit(dnds.dat)
    meth.frac = dnds.dat2[dnds.dat2[,Xcol] >= SEPARATOR,]
    un.meth.frac = dnds.dat2[dnds.dat2[,Xcol] < SEPARATOR,]
    mn.meth = mean(meth.frac[,Ycol])
    mn.un.meth = mean(un.meth.frac[,Ycol])
    print(paste("Mean Methylated =", mn.meth))
    print(paste("Mean Unmethylated =", mn.un.meth))
    print(paste("Percent Difference = ", (((mn.un.meth - mn.meth)/mn.meth)*100), "%", sep = ""))
    percentages = append(percentages, ((mn.un.meth - mn.meth)/mn.meth)*100)
    window_data = window_plot(1/8, Xcol, Ycol, dnds.dat2, Cutoff)
#     print("PRINTING WINDOW DATA:")
#     print(window_data)
    if (i < 2){
      plot(mn ~ x, type = "n", xlim = XLIM, ylim = mut.ylim, data = window_data, axes = F, xlab = "log2FoldChange", ylab = YLAB)
      axis(1)
      axis(2)
      loess_fit <- loess(mn ~ x, window_data, span = .9, se = T)
      lines(window_data$x, predict(loess_fit), lwd=LWD, col = line.color)
    }
    if (i >= 2){
      loess_fit <- loess(mn ~ x, window_data, span = .9, se = T)
      lines(window_data$x, predict(loess_fit), lwd=LWD, col = line.color)
    }
  }
print(paste("OVERALL MEAN DIFFERENCE FOR THIS SPECIES SET = ", mean(percentages)))
}

sub.species = 'Adigitifer'
for (sub.species in acro){
  dnds.dat = merge.dats(mdat, read.table(paste("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_", sub.species, '.txt', sep = ""), header = T), orthos, iso2seq, MUT.TYPE, FILT)
  print(sub.species)
  dnds.dat = na.omit(dnds.dat)
  print(dnds.dat[dnds.dat$dnds == max(dnds.dat$dnds),])
}
dnds.dat = dnds.dat[dnds.dat[,MUT.TYPE] != Inf,]


######### PLOT FOR ACROPORID COMPARISONS ###############
#dN
quartz()
LWD = 4
SEPARATOR = 0
par(mfrow = c(1, 3))
tot = sum(c(length(acro), length(non.acro), length(actin)))
COLOR.SET = rainbow(tot)
COLS = COLOR.SET[1:4]
MUT.TYPE = 'dN'
MUT.YLIM = c(0.006, 0.0245)
YLAB = "dN"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dS'
MUT.YLIM = c(.025, .14)
YLAB = "dS"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dNdS
MUT.TYPE = 'dnds'
MUT.YLIM = c(0.2, 0.47)
YLAB = "dN/dS"
plot.mut.stackF(acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
window_data = window_plot(1/8, 'log2FoldChange', MUT.TYPE, dnds.dat, Cutoff)

############### PLOT FOR OTHER CORALS ###############
#dN
###PLOT dN FOR NON ACROPORID SPECIES
quartz()
par(mfrow = c(1, 3))
COLS = COLOR.SET[6:9]
MUT.TYPE = 'dN'
MUT.YLIM = c(0.09, 0.155)
YLAB = "dN"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dS'
MUT.YLIM = c(1.25, 3.2)
YLAB = "dS"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dN/dS
MUT.TYPE = 'dnds'
MUT.YLIM = c(0.048, 0.1)
YLAB = "dN/dS"
plot.mut.stackF(non.acro, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

########## PLOT dS FOR ACTINIARIA #################
#dN
COLS = c('brown', 'tan')
quartz()
par(mfrow = c(1, 3))
MUT.TYPE = 'dN'
MUT.YLIM = c(0.18, 0.232)
YLAB = "dN"
plot.mut.stackF(actin, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dS'
MUT.YLIM = c(24.5, 30.5)
YLAB = "dS"
plot.mut.stackF(actin, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)
#dS
MUT.TYPE = 'dnds'
MUT.YLIM = c(.016,.021)
YLAB = "dN/dS"
plot.mut.stackF(actin, MUT.YLIM, MUT.TYPE, 'log2FoldChange', XLAB, YLAB, species, TRUE, 0, c(2, -2), COLS)

#PLOT THE LOESS LINES FOR PAIR-WISE MUTATION RATES FOR THE ACROPORID SPECIES ########

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

merge.cpg.dats = function(compare.species, SPECIES, mut.column){
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
  
  print(paste("Number of Genes before Outlier Removal = ", nrow(result)))
  result2 = result[result[,mut.column] < (mean(result[,mut.column]) + SDs*sd(result[,mut.column])), ]
  print(paste("Number of Genes after Outlier Removal = ", nrow(result)))
  head(result)
  return(result)
}

plot.cpg.mut = function(compare.species, speciesList){
  for (sp in speciesList){
    if (sp == compare.species){
      next
    }
    dat = merge.cpg.dats(compare.species, sp, MUT.TYPE)
    title <- paste(compare.species, 'vs', sp)
    plot.mut(dat, MUT.TYPE, 'cpgOE', MUT.TYPE, 'cpgOE', title, FALSE, 0.65)
  }
}
speciesList <- c('Amillepora', 'Apalmata', 'Ahyacinthus', 'Atenuis', 'Pastreoide', 'Pcarnosus', 'Pdamicorni', 'Spistillat', 'Apallida', 'Nvectensis')

quartz()
par(mfrow = c(1,4))
for (ANCHOR in speciesList){
  quartz()
  par(mfrow = c(3,5))
  plot.cpg.mut(ANCHOR, speciesList)
}

merge.cpg.dats()

################ PLOT dN or dS FOR ANY COMPARISON SPECIES AND A CHOSEN FAMILY ######################
quartz()
par(mfrow = c(1, 4))
MUT.TYPE = 'dN'
FAMILY <- acro
FAMILY <- non.acro
FAMILY <- actin
COMPARE.SP = 'Pdamicorni'
plot.cpg.mut(COMPARE.SP, acro)


#RUN FOR SET OF SPECIES
for (ANCHOR in non.acro){
  quartz()
  par(mfrow = c(1, 4))
  plot.cpg.mut(ANCHOR, acro)
}

########################################################################################




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








