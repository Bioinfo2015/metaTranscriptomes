##------------------------------------------------------------------------------
##------------ cgp_and_mutation_rates.R-----------------------------------------
##------------------------------------------------------------------------------
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
  x$sp = sp
  x$rep.strat = species[i,'rep.strat']
  dat = rbind(dat,x)
}
head(dat)
#===============================================================================

##------------------------------------------------------------------------------
##------------- FILTER BASED ON SIZE AND CPGO/E --------------------------------
##------------------------------------------------------------------------------
#set filtering parameters
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
dat2 = filter(dat, minlength, maxlength, maxOE, minOE, "All Species")

#Contruct histogram for all species
hist(dat2$cpgOE, breaks = 50, main = "All Species", prob = T, ylim = c(0, 1.4), xlab = "CpGoe", xlim = c(0,1.5))
#contruct histograms for brooders and spawners
# hist(dat2[dat2$rep.strat == 'brooder',]$cpgOE, breaks = 50, main = 'Brooders', prob = T, ylim = c(0, 1.4), xlab = "CpGoe", xlim = c(0,1.5), col = 'green')
# hist(dat2[dat2$rep.strat == 'spawn',]$cpgOE, breaks = 50, main = "Spawners", prob = T, ylim = c(0, 1.4), xlab = "CpGoe", xlim = c(0,1.5), col = 'blue')

##plot the histogram for each individual species
for (i in 1:length(species[,1])){
	sp = as.character(species[i,1])
	print(sp)
	st = species[i,2]
	sub = dat2[dat2$sp == sp,]
	hist(sub$cpgOE, breaks = 50, prob = T, main = paste(sp, st, nrow(sub), 'genes'))
}
#output the dataset
# write.table(cgm,'/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores7-27-14__1000LengthFilter.out',row.names=FALSE, quote = F)
##==============================================================================##
##==============================================================================##
##----------- LOOK AT RELATIONSHIP WITH dN and dS pairwise comparisons --------
#upload dnds data output data parsed from codeml
dnds = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/pair-wise_dNdS_digitifera.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_digitifera.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_digitifera11-25-14', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_digitiferaWanemone.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_Mcavernosa.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_siderea.txt', header = T)
head(dnds)

#upload the table of orthologs generated using get_reciprocal_orthos.py and mergeReciprocalOrthos.py 
orthos = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/Adigitifera_CDS_Ortholog_Table_e5_hp50_c3.txt', header = T)
colnames(orthos)[1] = "EST"
head(orthos)

#SET UP FUNCTION TO PLOT REGRESSIONS AND DO SOME STATS
require(scales)
plot_mut = function(mut.dat, TITLE){
  plot(dN~cpgOE, data = mut.dat, main = TITLE, xlim = rev(range(mut.dat$cpgOE)))
  lm1 = lm(dN~cpgOE, data = mut.dat)
  print(summary(lm1))
  abline(lm1, lwd = 2, col = 'red')
  spearmans.cor = cor(mut.dat$dN, mut.dat$cpgOE, method = "spearman")
  spear.test = cor.test(mut.dat$dN, mut.dat$cpgOE, method = "spearman")
#   print(spear.test)
}
# quartz()
# plot_mut(poc, 'acro')

#SET UP FUNCTION TO DIVIDE DATASETS INTO WINDOWS
#separates a dataset into quantiles based on CpGoe 
#and then plots the mean values for the column given
#for the argument 'column'
#quantile size is given by argument 'size'
library(plotrix)
quantile_plot = function(size, dat, Xcolumn, Ycolumn){
  windows = quantile(dat[,Xcolumn], probs = seq(0,1, size))
  #print(windows)
  mn = c()
  x = c()
  n = c()
  n.sig = c()
  strr = c()
  sig.ratio = c()
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat[,Xcolumn]>left,]
    sub = sub[sub[,Xcolumn]<right,]##subset out the percentile
    sub.n = nrow(sub)
    n = append(n, sub.n)
    mn = append(mn,mean(sub[,Ycolumn]))
    strr = append(strr, std.error(sub[,Ycolumn]))
    x = append(x, mean(sub[,Xcolumn]))
  }
  plot_dat = data.frame(mn,x,strr,n)
  return(plot_dat)
}
window_plot = function(size, slide, dat, column){
  lefts = seq(0,2, by = slide)
  rights = lefts + size
  #print(windows)
  mn = c()
  x = c()
  n = c()
  n.sig = c()
  strr = c()
  sig.ratio = c()
  for (i in 1:(length(lefts)-1)){
    left = lefts[i]
    right = rights[i]
    sub = dat[dat$cpgOE>left,]
    sub = sub[sub$cpgOE<right,]##subset out the percentile
    sub.n = nrow(sub)
    n = append(n, sub.n)
    mn = append(mn,mean(sub[,column]))
    strr = append(strr, std.error(sub[,column]))
    x = append(x, mean(sub$cpgOE))
  }
  plot_dat = data.frame(mn,x,strr,n)
  return(plot_dat)
}

#set up taxon sets
species
acroporidae = read.table("acroporidae.txt", col.names = c('sp'))
faviidae = read.table("faviidae.txt", col.names = c('sp'))
pocilloporidae = read.table("pocilloporidae.txt", col.names = c('sp'))
robust = read.table("robust.txt", col.names = c('sp'))
complex = read.table("complex.txt", col.names = c('sp'))
amillepora = read.table("amil.txt", header = T)


#RUNT THE FUNCTIONS ABOVE TO REFORMAT AND PLOT DATA FOR EACH SPECIES
#ALSO OUTPUTS EACH SPECIES AS ITS OWN DATAFRAME
#filter = the maximum mutation rate you will accept
#mut.type = the column heading of the data you want to plot 
#taxon = a list of species
# quartz()
reformat_mut_dat = function(taxon, filter, mut.type){
  spp = taxon
  print(spp)
  sub = dnds[dnds$species == spp,] ##get the substring becasue raxml uses onl first 10 characters
  head(sub)
  sub2 = merge(sub, orthos, by = 'EST')
  head(sub2)
  sub2$EST = sub2[,spp]
  head(sub2)
  sub3 = merge(sub2, dat2, by = 'EST')
  head(sub3)
  sub3 = sub3[sub3[,mut.type] < filter,]
  head(sub3)
  combined = rbind(combined, sub3)
  head(combined)
  print(paste(nrow(sub3), "Genes"))
  if (nrow(sub3) > 100){
    x = quantile_plot(1/12, sub3, 'cpgOE', mut.type)
    gene.count = sum(x$n)
    loess_fit <- loess(x$mn ~ x$x, span = 2, se = T)
    p = plotCI(x$x, x$mn, uiw = x$strr, err="y",xlab = expression("CpG"["O/E"]), pch=19, cex = .5, lwd=2, sfrac = .005, main=paste(mut.type, spp, gene.count, 'genes'), axes = T, cex.lab = 1.1, add = F, col = 'black')
    lines(x$x, predict(loess_fit),col="red",lwd=1)
  }
  assign(paste(spp, '.mut', sep = ''), sub3)
#   return(combined)
}
reformat_mut_dat("Amillepora", 1, 'dN')
x=reformat_mut_dat("Mauretenra", 1, 'dN')

quartz()
par(mfrow=c(1,1))
FILT = 1
MUT = 'dN'
combined = reformat_mut_dat(species, FILT, MUT)
acro = reformat_mut_dat(acroporidae, FILT, MUT)
rob = reformat_mut_dat(robust, FILT, MUT)
comp = reformat_mut_dat(complex, FILT, MUT)
fav = reformat_mut_dat(faviidae, FILT, MUT)
poc = reformat_mut_dat(pocilloporidae, FILT, MUT)
amil = reformat_mut_dat(amillepora, FILT, MUT)
amil = reformat_mut_dat(amillepora, FILT, MUT)

quartz()
for(i in species$sp){
  spp = as.character(i)
  reformat_mut_dat(spp, FILT, MUT)
}



dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_Nventensis.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_digitiferaWanemone.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_Mcavernosa.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_Mauretenra.txt', header = T)
dnds = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/pair-wise_dNdS_siderea.txt', header = T)

##----- compare dN and dS with expression levels from Carly's data ----------
tdat = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/cpg_OE/expressionData.txt", header = T)
iso2seq = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/mgi_GO_slim/amillepora_transcriptome_july2014/amil_seq2iso.tab", col.names = c('EST', 'isogroup'))
amil = acro[acro$species == 'Amillepora',]
acro$species[acro$species == 'Amillepora']
head(amil)
head(iso2seq)
x = merge(amil, iso2seq, by = 'EST')
head(x)
x$EST = x$isogroup
y = merge(x, tdat, by = 'EST')
nrow(y)
head(y)
plot(dS~grandmeans, data = y, ylim = c(0,5))
lm2 = lm(dS~grandmeans, data = y)
summary(lm2)
abline(lm2, col = 'red')

mut.type = 'dS'
spp = 'A.millepora'
x = quantile_plot(1/12, y, 'grandmeans', mut.type)
gene.count = sum(x$n)
loess_fit <- loess(x$mn ~ x$x, span = 2, se = T)
p = plotCI(x$x, x$mn, uiw = x$strr, err="y",xlab = expression("expression"), pch=19, cex = .5, lwd=2, sfrac = .005, main=paste(mut.type, spp, gene.count, 'genes'), axes = T, cex.lab = 1.1, add = F, col = 'black', xlim = rev(range(x$x)))
lines(x$x, predict(loess_fit),col="red",lwd=1)

##----------- COMPARE RESULTS FOR CPGOE BETWEEN TRANSCRIPTOME ONLY AND GENE BODIES FROM GENOME --------
##the idea behind this section is that we are switching to just transcriptome.
##to set this up, first run the sections above for gene bodies then set cgmGeneBodies = cgm, then rerun for just transcriptome
# cgmGeneBodies = cgm
# compare = merge(cgmGeneBodies, cgm, by = 'EST')
# head(compare)
# plot(cpgOE.x~cpgOE.y, data = compare)
# lm.test = lm(cpgOE.x~cpgOE.y, data = compare)
# summary(lm.test)
#exporting the cpg scores:
# write.table(cgm,'/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/cpgData_Routput.txt',row.names=FALSE, quote = F)
##--------------------- MIXTURE MODELING ---------------------------------------
##use this section to manipulate the model parameters and find the best ones the best parameters are hard-coded in the section above so that figures can be made with set colors.
library(mixtools)
emmix = function(cg,species,means,sigma,lambda,k){
  x = cg$cpgOE
  cpgNull <- normalmixEM(x, mu = means, sigma = sigma, lambda = lambda, k=k,epsilon=0.01, arbvar=T)
  summary(cpgNull)
  plot(cpgNull, which = 2, breaks = 30, density = TRUE, cex.axis = 1.4, cex.lab = 1.5, cex.main = 1.5, main2 = species, xlab2 = "CpGo/e") 
  return(cpgNull)
}
par(mfrow=c(1,1))
##NULLs
mil2mix = emmix(cgm,"millepora",NULL,NULL,NULL,2)
#dig2mix = emmix(cgd,"digitifera",NULL,NULL,NULL,2)
model_df = function(lambda,mean,sigma){
  df = data.frame(cbind(lambda,mean,sigma))
  rownames(df) = paste("comp",rownames(df),sep="")
  return(df)
}##assigns the mixture model to a variable
mmix2 = model_df(mil2mix$lambda,mil2mix$mu,mil2mix$sigma)
##============================================================================
##------------ TESTING FOR NUMBER OF COMPONENTS --------------
par(mfrow = c(1,1))
x.norm = pnorm(cgm$cpgOE, mean = 0.704, sd = 0.2394679)
x = sample(cgm$cpgOE, 5000, replace = F)
shapiro.test(x)
#reject the null hypothesis that the distribution is normal
require(mclust)
x = cgm$cpgOE
clust = Mclust(x, G = c(1:3), modelNames = c("V"))
bic = mclustBIC(x, G = c(1:5), modelNames = c("V"))
plot(bic, G = c(1:5), modelNames = c("V"))