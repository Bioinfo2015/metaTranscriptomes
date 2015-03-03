#This R script is for performing likelihood ratio tests on outputs from paml
#It is specifically designed to compare a Null model where dN/dS (W) is equal across branches
#against and alternative model where W is greater in the Acroporid lineage
#-----------------------------------------------------------------------
#-------------EXAMPLE OF LIKELIHOOD RATIO TEST--------------------------
#-----------------------------------------------------------------------
#INFORMATION TAKEN FROM WHITLOCK AND SCHLUTER (ANALYSIS OF BIOLOGICAL DATA 2009)
#PICTURE TWO MODELS, A NULL MODEL (Mo) REPRESESENTING THE NULL HYPOTHESIS (Ho) 
#AND AN ALTERNATIVE MODEL (Ma) REPRESENTING THE ALTERNATIVE HYPOTHESIS (Ha)
#EACH HAS AN ESTIMATED LIKELIHOOD (from separate codeml runs)
#WHICH DESCRIBES THE PROBABILITY OF GETTING OUR DATA GIVEN THE MODEL PARAMETERS
#START WITH THESE VALUES
L.Mo = 1e-10
L.Ma = 1e-9
#CODEML GIVES LIKELIHOOD AS LOG LIKELIHOODS, SO IT WILL REALLY BE LIKE THIS
lnL.Mo = log(L.Mo)
lnL.Ma = log(L.Ma)
#CODEML WILL ALSO OUTPUT THE NUMBER OF FREE PARAMETERS FOR EACH MODEL
#IN THIS CODEML OUTPUT THE VALUES FOR LIKELIHOOD AND FREE PARAMETERS
#LOOKS LIKE THIS: lnL(ntime: 25  np: 28):  -4843.605533
#WHERE np: 28 GIVES FREE PARAMETERS AND -4843.605533 IS THE LOG LIKELIHOOD
np.Mo = 28
np.Ma = 29
#NOW WE HAVE ALL THE DATA WE NEED TO PERFORM A LOG LIKELIHOOD TEST
#GET THE LIKELIHOOD RATIO TEST STATISTIC (G)
G = 2*(lnL.Ma - lnL.Mo) 
#(in terms of likelihood it would look like this: G = 2*log(likelihood.Ma/likelihood.Mo) )
#WHEN THE NULL HYPOTHESIS IS TRUE, G HAS AN APPROXIMATE CHI-SQUARE DISTRIBUTION
#SO WE CAN USE THE CHI-SQUARE DISTRIBUTION WITH THE RIGHT NUMBER OF DEGREE OF FREEDOM (df)
#TO CALUCLATE A P VALUE FROM THE G STATISTIC
#DEGREES OF FREEDOM IS EQUL TO THE DIFFERENCE BETWEEN THE TWO HYPOTHESES IN THE NUMBER OF PARAMTERS
#THAT WERE ESTIMATED FROM THE DATA (this is our np value for each model)
df = np.Ma - np.Mo
#NOW WE CAN EXTRACT A P VALUE FROM THE CHI-SQUARE DISTRIBUTION
#THIS REPRESENTS
p.value = pchisq(G, df, ncp = 0, lower.tail = F)
p.value

#TO SEE HOW THIS WORKS WE CAN PLOT THE CHI-SQUARE PROBABILITY DENSITY CURVE AND SEE WHERE OUR G IS
x <- seq(-.1,10,by=.001)
y <- dchisq(x, df=1)
plot(x,y, pch = '.')
abline(v = G)
dchisq(G, df = 1)

x = seq(from = 0, to = 8, by = 0.001)
curve (dchisq(x, df = 5))
#use help(Chisquare) for more information
#FOR STARTING LIKELIHOODS L.Mo = 1e-10 and L.Ma = 1e-9 WE GET P-VALUE = 0.0024
#BASED ON THIS WE REJECT THE NULL MODEL. THE ALTERNATIVE MODEL PROVIDES SIGNIFICANTLY BETTER FIT TO OUR DATA.
#HERE IS A FUNCTION TO PEFORM THIS
lrt = function(lnL.Ma,  lnL.Mo, np.Mo, np.Ma){
  G = 2*(lnL.Ma -  lnL.Mo)
  df = np.Ma - np.Mo
  p.value = pchisq(G, d.f., ncp = 0, lower.tail = F)
  print(p.value)
}
#=============================================================================

#-----------------------------------------------------------------------
#------------- RUNNING LRTs FOR EACH ORTHOLOG --------------------------
#-----------------------------------------------------------------------
#NOW START WITH A SET OF NULL AND ALTERNATIVE MODELS FOR A SET OF ORTHOLOGS
#OUR QUESTION IS WHICH ORTHOLOGS ARE BETTER DESCRIBED WITH OUR ALTERNATIVE MODEL

#READ IN THE SET OF MODEL STATISTICS
setwd("/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis")
dat = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/data_analysis/likelihoods.txt', header = T)
head(dat)
?apply

#SET UP FUNCTION TO RETURN THE P VALUE FOR THE LRT BETWEEN TWO CODEML MODELS
lrt = function(vector){
  lnL.Mo = vector[1]
  lnL.Ma = vector[2]
  np.Mo = vector[3]
  np.Ma = vector[4]
  G = 2*(lnL.Ma -  lnL.Mo)
  df = np.Ma - np.Mo
  p.value = pchisq(G, df, lower.tail = F)
}
dat$p.values = apply(dat[,2:5], 1, lrt)
head(dat)

#UPLOAD MBD-SEQ COUNT DATA SCORES
mbdat = read.table('/Users/grovesdixon/Documents/lab_files/metaTranscriptomes/pull_down/mbdScoreIsogroups.txt', header = T)
EST = mbdat$EST
mbd.score = mbdat$mbd.score
cpg = mbdat$cpgOE
mbdat = data.frame(EST, mbd.score, cpg)
head(mbdat)

#UPLOAD DATASETS NEEDED TO TRANSITION BETWEEN GENE NAMES
#UPLOAD MILLEPOR-DIGITIFERA ORHTOLOG TABLE
odat = read.table('/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/Adigitifera-REF_Amillepora_Protein_Reciprocal_Blast_Hits_12-20-14.txt', header = T)
colnames(odat) = c("dig.names", "EST")
head(odat)
#UPLOAD THE MILLEPORA SEQ2ISO TABLE TO GO FROM CONTIGS TO ISGROUPS
seq2iso = read.table('/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/amillepora_transcriptome_july2014/amil_seq2iso.tab', header = F)
colnames(seq2iso) = c('contig', 'EST')
head(seq2iso)
#NOW DO A LOT OF MERGING TO GET A CONCENSUS DATAFRAME USING 'EST' AS THE MAIN HANDLE FOR GENE NAMES
change.to.contig.names = merge(mbdat, seq2iso, by = 'EST')
change.to.contig.names$EST = change.to.contig.names$contig
head(change.to.contig.names)
change.to.dig.names = merge(change.to.contig.names, odat, by = 'EST')
change.to.dig.names$EST = change.to.dig.names$dig.names
head(change.to.dig.names)
#MERGE WITH THE LIKELIHOOD RATIO TEST DATAFRAME FROM ABOVE
head(dat)
colnames(dat) = c('EST', 'nl', 'al','nn','an','p.values')
dat3 = merge(change.to.dig.names, dat, by = 'EST')
head(dat3)
dat4 = data.frame(dat3$EST, dat3$mbd.score, dat3$cpg, dat3$p.values)
colnames(dat4) = c("EST", "mbd.score", "cpgoe", "p.values")
dat4$adj.p = p.adjust(dat4$p.values, method = 'BH')
head(dat4)
#CLEAN UP THE MESS WE MADE GETTING THERE
rm(change.to.contig.names)
rm(change.to.dig.names)
rm(odat)
rm(seq2iso)
rm(dat3)


#CORRELATION BETWEEN POSITIVE SELECTION IN ACROPORIDS AND MBD-SCORE
head(dat4)
hist(dat4$mbd.score, breaks = 50)
dat4$adj.p = p.adjust(dat4$p.values, method = 'BH')
plot(-log(adj.p)~mbd.score, data = dat4)
lm1 = lm(-log(p.values)~mbd.score, data = dat4)
summary(lm1)
abline(lm1)


require(plotrix)
size = .2
dat = dat4
Xcolumn
quantile_plot = function(size, dat, Xcolumn, Ycolumn, sig.threshold){
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
    sub1 = dat[dat[,Xcolumn]>left,]
    sub= sub1[sub1[,Xcolumn]<right,]##subset out the percentile
    sub.n = nrow(sub)
    sig = sub[sub[,Ycolumn] < sig.threshold,]
    sig.count = nrow(sig)
    n.sig = append(n.sig, sig.count)
    n = append(n, sub.n)
    mn = append(mn,mean(sub[,Ycolumn]))
    strr = append(strr, std.error(sub[,Ycolumn]))
    x = append(x, mean(sub[,Xcolumn]))
  }
  plot_dat = data.frame(mn,x,strr,n, n.sig)
  return(plot_dat)
}

qdat = quantile_plot(1/20, dat4, 'mbd.score', 'adj.p', 0.001)
qdat = quantile_plot(1/20, dat4, 'cpgoe', 'adj.p', 0.001)
plot(n.sig~x, data = qdat)
loess_fit <- loess(n.sig ~ x, qdat, span = 1.5, se = T)
lines(qdat$x, predict(loess_fit),col="red",lwd=1)






