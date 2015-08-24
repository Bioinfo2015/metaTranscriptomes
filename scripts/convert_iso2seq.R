#convert_iso2seq.R
#This script takes in a list of isogroups from the A.millepora transcptome and returns the corresponding list of contigs based on the seq2iso table

#upload seq2iso table
iso2seq = read.table("/Users/grovesdixon/Documents/lab_files/Amillepora_transcriptome/amillepora_transcriptome_july2014/amil_seq2iso.tab", col.names = c('EST', 'isogroup'))
head(iso2seq)

#get file info and upload it
path = readline(prompt="Enter path to the file: ")
fileName = readline(prompt='Enter the file name: ')
isos = read.table(paste(path, fileName, sep = "/"))
colnames(isos) = c('isogroup')
head(isos)

#merge the tables to get contigs
dat = merge(isos, iso2seq, by = 'isogroup')
head(dat)

#output the file
outName = readline(prompt("What output name? "))
write.table(dat$EST, paste(path, outName, sep = "/"), row.names = F, quote = F)
