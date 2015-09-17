#convert_iso2seq.R
#This script takes in a list of isogroups from the A.millepora transcptome and returns the corresponding list of contigs based on the seq2iso table

#upload seq2iso table
setwd("/Users/grovesdixon/git_Repositories/metaTranscriptomes/working_directory")
iso2seq = read.table("amil_seq2iso.tab", col.names = c('contig', 'isogroup'))
head(iso2seq)

#get file info and upload it
isos = read.table('ribosomalIsogroups.txt')
colnames(isos) = c('isogroup')
head(isos)

#merge the tables to get contigs
dat = merge(isos, iso2seq, by = 'isogroup')
head(dat)

#output the file
write.table(dat$contig, 'ribosomalContigs.txt', row.names = F, quote = F)
