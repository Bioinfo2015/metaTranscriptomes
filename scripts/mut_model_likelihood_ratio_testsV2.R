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

x = c(np.Ma, lnL.Ma, np.Mo, lnL.Mo)
lrt(x)
#=============================================================================

########################################################################################
############################ RUNNING LRTs FOR EACH ORTHOLOG ############################
########################################################################################

#NOW START WITH A SET OF NULL AND ALTERNATIVE MODELS FOR A SET OF ORTHOLOGS
#OUR QUESTION IS WHICH ORTHOLOGS ARE BETTER DESCRIBED WITH OUR ALTERNATIVE MODEL

#READ IN THE SET OF MODEL STATISTICS
setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files")
altdat = read.table('altLikelihoods.txt', header = T)
nulldat = read.table('nullLikelihoods.txt', header = T)
head(altdat)
head(nulldat)
ldat = merge(altdat, nulldat, by = 'EST')
colnames(ldat) = c('EST', 'np.a', 'l.a', 'np.o', 'l.o')
nrow(ldat)
head(ldat)

#SET UP FUNCTION TO RETURN THE P VALUE FOR THE LRT BETWEEN TWO CODEML MODELS
ldat$G = 2*(ldat$l.a - ldat$l.o)
ldat$df = ldat$np.a - ldat$np.o
ldat$p.value = pchisq(ldat$G, ldat$df, lower.tail = F)
ldat$adj.p = p.adjust(ldat$p.value, method = 'BH')


#UPLOAD MBD-SEQ COUNT DATA SCORES
mdat = read.table('methylDESeq_v1.txt', header = T)
mdat$isogroup = rownames(mdat)

#UPLOAD DATASETS NEEDED TO TRANSITION BETWEEN GENE NAMES
#UPLOAD MILLEPOR-DIGITIFERA ORHTOLOG TABLE
#UPLOAD THE ORTHOLOG DATA
orthos = read.table('/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/data_files/Adigitifera_CDS_Ortholog_Table_e5_hp50_c2.txt', header = T)##best
#UPLOAD THE SEQ2ISO DATA FOR MILLEPORA
iso2seq = read.table("/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/amillepora_transcriptome_july2014/amil_seq2iso.tab", col.names = c('EST', 'isogroup'))
colnames(iso2seq) = c('Amillepora', 'isogroup')

colnames(orthos)[1] = 'EST'

#WE NEED THE MBD SEQ DATA COMBINED WITH THE MUTATION MODEL DATA
#THE SUBSTITUTION DATASET IS IN TERMS OF A.DIGITIFERA SEQUENCES
#BUT THE MBD-SEQ DATA ARE IN TERMS OF A.MILLEPORA ISOGROUPS
#SO WE NEED TO USE THE ORTHOLOG TABLE TO CONVERT THE NAMES AROUND
#BUT THE ORTHOLOG TABLE HAS A.MILLEPORA CONTIG NAMES INSTEAD OF ISOGROUP NAMES
#SO SEVERAL DATASETS MUST BE MERGED TO GET EVERYTHING IN THE RIGHT PLACE
#first we need the millepora ortholog data in 'isogroups' instead of 'contigs'
#so merge the ortholog data with the iso2seq data
head(iso2seq)
head(orthos)
m1 = merge(iso2seq, orthos, by = 'Amillepora')
head(m1)
#we only need the contigs, isogroups, and Adigitifera seqs
m1 = m1[,1:3]
head(m1)
#now pull the likelihood data in
m2 = merge(m1, ldat, by = 'EST')
head(m2)
#now pull in the MBD-seq data
cdat = merge(m2, mdat, by = 'isogroup')
head(cdat)

#CLEAN UP THE MESS WE MADE GETTING THERE
rm(m1)
rm(m2)
rm(iso2seq)
rm(orthos)

#WHAT ARE THE GENES THAT ARE UNDER SELECTION?
CUT = 0.05
nrow(cdat[cdat$adj.p < CUT,])
cdat[cdat$adj.p < CUT,]
write.table(cdat[cdat$adj.p < CUT, 'EST'], 'isogroupsUnderSelection.txt', sep = '\t', quote = F, row.names = F)
nrow(cdat[cdat$p.value < CUT,])

#CORRELATION BETWEEN POSITIVE SELECTION IN ACROPORIDS AND MBD-SCORE
plot(-log(p.value, 10)~log2FoldChange, data = cdat)
lm1 = lm(-log(p.value, 10)~log2FoldChange, data = cdat)
summary(lm1)
abline(lm1)

require(plotrix)
size = .2
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
qdat = quantile_plot(1/20, cdat, 'log2FoldChange', 'p.value', 0.1)
plot(n.sig~x, data = qdat)
loess_fit <- loess(n.sig ~ x, qdat, span = 1.5, se = T)
lines(qdat$x, predict(loess_fit),col="red",lwd=1)
#SO APPARENTLY NO CORRELATION WHATSOEVER WITH MBD-SCORE

################# NOW SAVE THE DATAFRAME TO DO SOME GO ENRICHMENT ####################
#MISHA'S SCRIPTS USE ISOGROUP100 INSTEAD OF ISOGROUP=100, SO CHANGE THAT HERE
new.est = c()
for (i in cdat$isogroup){
  y = paste(strsplit(i, '=')[[1]][1], strsplit(i, '=')[[1]][2], sep = '')
  new.est = append(new.est, y)
}
out = data.frame(new.est, -log(cdat$p.value, 10))
colnames(out) = c('gene', 'logp')
head(out)
nrow(out)
#write out for GO
write.table(out, '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/go/Acroporid_evolution_pvalues2-21-15.txt', quote = F, row.names = F, sep = ",")
#write out for KOGG
write.table(out, '/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/kog/Acroporid_evolution_pvalues2-21-15.txt', quote = F, row.names = F, sep = ",")


