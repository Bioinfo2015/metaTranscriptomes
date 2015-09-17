setwd("/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory")
dat = read.table('pair-wise_dNdS_Amillepora.txt', header = T)
head(dat)


#
XLIM = c(0, 0.1)
spp = 'Adigitifer'
sub = dat[dat$species == spp,]
hist(sub$dNdS, xlim = XLIM)
hist(sub$dN, xlim = XLIM)
plot(density(sub$dN), xlim = XLIM)
x = sub[sub$dNdS > 1,]


#function to remove genes with dS values equal to zero since these provide meaningless dN/dS values
#also filters genes with dN less than the given value for dNcut
filter = function(dat, dNcut){
  #remove genes with no synonymous subs
  f = dat[dat$dS > 0,]
  #remove genes with dN values lower than a given threshold
  f = f[f$dN < dNcut,]
  return(dat)
}

#fuction get.candidates
#return a set of candidate genes under positive selection
#looks for genes that have dNdS above one for both pairwise comparisons
get.candidates = function(dat, sppa, sppb){
  #first subset for comparisons against each other species
  a = dat[dat$species==sppa,]
  b = dat[dat$species==sppb,]
  #remove any genes with dNdS above 2, because these likely represent bad ortholog calls
  print("This many genes have dNdS above 2 and are likely spurious orhtolog calls:")
  print(nrow(a[a$dNdS >= 2,]))
  print(nrow(b[b$dNdS >= 2,]))
  a = a[a$dNdS < 2,]
  b = b[b$dNdS < 2,]
  pa = a[a$dNdS > 1,]
  pb = b[b$dNdS > 1,]
  o = merge(pa, pb, by = 'contig')
  print("Number of Candidate Genes Under Positive Selection:")
  print(nrow(o))
  return(o)
}
p.carb = get.candidates(cdat, 'flav', 'macro')
p.flav = get.candidates(fdat, 'carb', 'macro')
p.macro = get.candidates(mdat, 'flav', 'carb')













