#R script LRT_for_sites_models.R
#read and analyze data from all five sites models in CODEML (model = 0, NSsites = 0, 1, 2, 7, 8)
#the different sites models are refered to by the value given for NSsites in the codeml control file


#read in the data
setwd("/Users/grovesdixon/Documents/lab_files/coding4people/moises/data/")
dat = read.table("sites_models_results_Rinput.txt", header = T)
tot = nrow(dat)
head(dat)

#set up some global variables
INITIAL.CUT = 0.05 #the p value cutoff for keeping genes after first test
SIG.CUT = 0.001    #the adjusted p values threshold to use for significance

#set up lrt function that returns p values from two vectors of likelihoods
lrt = function(la, lo, df){
  G = 2*(la - lo)
  p.values = pchisq(G, df, ncp = 0, lower.tail = F)
}

#first run lrt between model 0 and model 1
#filter the dataset to only include genes that are somewhat significant
p.values = lrt(dat$l1, dat$l0, 1)
dat$test1.p = p.values

#set up second reduced dataset only including genes that were significant for first test
dat2 = dat[dat$test1.p < INITIAL.CUT,]
nrow(dat2)

#run lrt between model 1 and 2
p.values = lrt(dat2$l2, dat2$l1, 2)
dat2$test2.adjp = p.adjust(p.values, method = "BH")
sig = dat2[dat2$test2.adjp < SIG.CUT,]
nrow(sig)
nrow(sig)/tot
head(sig)

#run lrt between model 7 and 8
dat2$test3.adjp = p.adjust(lrt(dat2$l8, dat2$l7, 2), method = "BH")
sig3 = dat2[dat2$test3.adjp < SIG.CUT,]
nrow(sig3)

#check the overlap
sig.both = merge(sig, sig3, by = 'gene')
nrow(sig.both)

