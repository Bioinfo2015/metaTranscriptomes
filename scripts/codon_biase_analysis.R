#codon_biase_analysis.R
#Groves Dixon
#9/1/15
#This script tests hypotheses of how gene body methylation should influence 
#codon bias data output from codonW

#1. Arginine codons lacking CpGs should be preferred
#2. 


#set working directory
setwd('/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory')


#upload and sort the data
#hilo data
hdat = read.table('hilo_data.txt', col.names = c('codon', 'aa', 'hi.use', 'hi.count', 'low.use', 'low.count', 'optimal.codon', 'cg.codon', 'gc.codon', 'total.hi', 'total.low', 'species'))
hdat = hdat[order(hdat$species, hdat$aa), ]
head(hdat)

#cai data
cdat = read.table('cai_data.txt', col.names = c('codon', 'aa', 'xi', 'wi', 'cg.codon', 'gc.codon', 'species'))
cdat = cdat[order(cdat$species, cdat$aa),]
head(cdat)

#set up species list
spp = levels(hdat$species)

#test hypothesis 1--Non CpG Arganines are preferred
optimal = c()
for (sp in spp){
	sub = hdat[hdat$species == sp,]
	ar = sub[sub$aa == 'Arg',]
	op = ar[ar$cg.codon == 0,]
	optimal = append(optimal, sum(op$optimal.codon))
}
res = data.frame(spp, optimal)







cgdat = hdat[hdat$cg.codon == 1,]
hi = cgdat$hi.use > cgdat$low.use
hi[hi == 'TRUE'] <- 1
hi[hi == 'FALSE'] <- 0
cgdat$hi = hi
wrong = cgdat[cgdat$optimal.codon == 1,]
right = cgdat[cgdat$optimal.codon == 0,]

#seems true for all except Apallida and Pdamicornis, where something weird is going on



sub = hdat[hdat$species == 'Adigitifera',]
sub = hdat[hdat$aa == 'Arg',]


