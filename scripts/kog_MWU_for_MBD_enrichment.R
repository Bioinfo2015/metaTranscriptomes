setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/kog")

# change these to your file names
inname = "MBD_foldChangeForKOG_V2_5-14-15.txt"
gene2kog="amil_defog_iso2kogClass.tab"

# change to Alternative="l" if the measure of interest gets more interesting as it gets smaller, i.e., p-value (the example brown.csv has -log() transformed p-values so in that case "g" is correct)
Alternative="t"

# mark everything (command-A) and execute!
#---------------------
mwut=function(gos) {
	terms=unique(gos$term)
	gos$seq=as.character(gos$seq)
	nrg=gos[!duplicated(gos$seq),2]
	names(nrg)=gos[!duplicated(gos$seq),1]
#	nrg=nrg+rnorm(nrg,sd=0.01) # to be able to do exact wilcox test
	rnk=rank(nrg)
	names(rnk)=names(nrg)
	pvals=c();drs=c();nams=c();levs=c();nseqs=c()
	for (t in terms){
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		sgo.yes=got$seq
		n1=length(sgo.yes)
		sgo.no=ngot$seq
		n2=length(sgo.no)
		if (n2 < n1) {
			print(paste("skipping",t,"nseqs =",n1))
			next
		}
		wi=wilcox.test(nrg[sgo.yes],nrg[sgo.no],alternative=Alternative)	# removed correct=FALSE 
		r1=sum(rnk[sgo.yes])/n1
		r0=sum(rnk[sgo.no])/n2
		dr=r1-r0
		drs=append(drs,round(dr,0))
		levs=append(levs,got$lev[1])
#		nams=append(nams,as.character(got$term[1]))
		pvals=append(pvals,wi$p.value)
		nseqs=append(nseqs,n1)	
	}
	res=data.frame(cbind(nseqs,"delta.rank"=drs,"pval"=pvals))
	res=cbind("term"=as.character(terms),res)
	res$pval=as.numeric(as.character(res$pval))
	res$delta.rank=as.numeric(as.character(res$delta.rank))
#	res$level=as.numeric(as.character(res$level))
	res$nseqs=as.numeric(as.character(res$nseqs))
	res=res[order(res$pval),]
	res$padj=p.adjust(res$pval,method="BH")
	return(res)
}

ft=function(gos) {
gos=annotated
	terms=unique(gos$term)
	gos$seq=as.character(gos$seq)
	pft=c();nam=c();lev=c();nseqs=c()
	for (t in terms) {
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		go.sig=sum(got$value)
		go.ns=length(got[,1])-go.sig
		ngo.sig=sum(ngot$value)
		ngo.ns=length(ngot[,1])-ngo.sig
		sig=c(go.sig,ngo.sig) # number of significant genes belonging and not belonging to the tested GO category
		ns=c(go.ns,ngo.ns) # number of not-significant genes belonging and not belonging to the tested GO category
		mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("go","notgo"),sig=c("go","notgo")))
		ff=fisher.test(mm,alternative="greater")
		pft=append(pft,ff$p.value)
		nam=append(nam,as.character(got$name[1]))
		lev=append(lev,got$lev[1])
		nseqs=append(nseqs,length(got[,1]))
	}
	res=data.frame(cbind("term"=as.character(terms),nseqs,"pval"=pft))
	res$pval=as.numeric(as.character(res$pval))
	res$nseqs=as.numeric(as.character(res$nseqs))
	res=res[order(res$pval),]
	res$padj=p.adjust(res$pval,method="BH")
	return(res)
}

#---------------------
rsq=read.csv(inname)
names(rsq)=c("seq","value")
kogs=read.table(gene2kog,sep="\t")
annotated=rsq[rsq$seq %in% kogs$V1,]
kogrows=match(annotated$seq,kogs$V1)
kogs=kogs[kogrows,]
annotated$term=kogs$V2
#head(annotated)

mwut.t=TRUE
if (length(levels(as.factor(annotated$value)))==2) {
	print("Binary classification detected; will perform Fisher's test");
	mwut.t=F
	rr=ft(annotated)
} else {
	print("Continuous measure of interest: will perform MWU test");		
	rr=mwut(annotated)
}
rr

write.table(rr, "KOG_Results_MBD-seq_Enrichment_5-15-15.txt", sep="\t", quote=F)



#--------------------- RUN FOR ALL SPECIES ----------------------
spp.list = c("Adigitifer", "Ahyacinthu", "Apalmata", "Atenuis", "Pastreoide", "Pcarnosus", "Pdamicorni", "Spistillat")
for (i in spp.list){
  inname = paste(i, "-Amillepora_dnds.txt", sep = "")
  rsq=read.csv(inname)
  names(rsq)=c("seq","value")
  kogs=read.table(gene2kog,sep="\t")
  annotated=rsq[rsq$seq %in% kogs$V1,]
  kogrows=match(annotated$seq,kogs$V1)
  kogs=kogs[kogrows,]
  annotated$term=kogs$V2
  #head(annotated)
  
  mwut.t=TRUE
  if (length(levels(as.factor(annotated$value)))==2) {
    print("Binary classification detected; will perform Fisher's test");
    mwut.t=F
    rr=ft(annotated)
  } else {
    print("Continuous measure of interest: will perform MWU test");		
    rr=mwut(annotated)
  }
  print(paste("Results for Comparison between", i, "and A. millepora:"))
  print(rr)
  
  write.table(rr, paste("KOG_results_Continuous_dNds_", i, "V_Amillepora_5-8-15.txt", sep = ""), sep="\t", quote=F)
}
