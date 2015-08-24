setwd("/Users/grovesdixon/Documents/lab_files/projects/metaTranscriptomes/go")


# Edit these as needed, then highlight everything and execute (takes ~5 min for MF and BP):
input = "MBD-scores_for_GO_MWU.csv" ## use this for go enrichment by MBD fold changes

goAnnotations="amil_defog_iso2go.tab"
goDatabase="go.obo.txt"
goDivision="MF"####MF is most interesting here
extraOptions="alternative=t"  #alternative=g/l/t
absValue=-log(0.1, 10) # we want to count genes with this or better absolute measure value within each displayed GO category. This does not affect statistics.

level1=1e-6 # p-value cutoff for the GO categories to plot as a tree; the worst ones will be printed in small italic font. Specify cutoff=1 to summarize all the genes at or exceeding absValue. 
level2=1e-8# # p-value cutoff for printing in regular font.
level3=1e-10 # p-value cutoff for printing in large bold font.
adjusted=T # replace with F to plot un-adjusted p-values.
txtsize=1  # decrease this one to squeeze more GO descriptions on the same panel.
font.family="sans" #"serif"

################################################################
# generating MWU and dissim files; only need to do this once
system(paste("perl ./gomwu.pl", goDatabase, goAnnotations, input, goDivision, extraOptions))
################################################################

require(ape)
in.mwu=paste("MWU", goDivision, input,sep="_")
in.dissim=paste("dissim", goDivision, goAnnotations,sep="_")

cutoff=-log(level1,10)
pv=read.table(in.mwu,header=T)
row.names(pv)=pv$term
in.raw=paste(goDivision,input,sep="_")
rsq=read.table(in.raw,sep="\t",header=T)
rsq$term=as.factor(rsq$term)

if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
heat=data.frame(cbind("pval"=pvals)) 
row.names(heat)=pv$term
heat$pval=-log(heat$pval+1e-15,10)
heat$direction=0
heat$direction[pv$delta.rank>0]=1
if (cutoff>0) { 
	goods=subset(heat,pval>=cutoff) 
} else {
	goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
	goods=heat[row.names(heat) %in% goods.names,]
}

colors=c("dodgerblue2","firebrick1","skyblue","lightcoral")
if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
	colors=c("black","black","grey50","grey50")
}
goods.names=row.names(goods)

# reading and subsetting dissimilarity matrix
diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
row.names(diss)=names(diss)
diss.goods=diss[goods.names,goods.names]

# how many genes out of what we started with we account for with our best categories?
good.len=c();good.genes=c()
for (g in goods.names) {
	sel=rsq[rsq$term==g,]	
	pass=abs(sel$value)>=absValue
	sel=sel[pass,]
	good.genes=append(good.genes,as.character(sel$seq))
	good.len=append(good.len,nrow(sel))
}
ngenes=length(unique(good.genes))

#hist(rsq$value)
totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")

# clustering terms better than cutoff
GO.categories=as.dist(diss.goods)
cl.goods=hclust(GO.categories,method="average")
labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
goods=goods[labs,]

quartz()
plot(as.phylo(cl.goods),show.tip.label=FALSE)

step=100
left=1
top=step*(2+length(labs))
quartz()
plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
ii=1
goods$color=1
goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
for (i in length(labs):1) {
	ypos=top-step*ii
	ii=ii+1
	if (goods$pval[i]> -log(level3,10)) { 
		text(left,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family) 
	} else {
		if (goods$pval[i]>-log(level2,10)) { 
			text(left,ypos,labs[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
		} else {
#			if (goods$pval[i]>cutoff) { 
#				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
	#		} else { 
		text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family) 
		#}
		}
	}
}

quartz()

plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
text(left,top,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
text(left,top-step,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
text(left,top-step*2,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)

print(paste("GO terms dispayed: ",length(goods.names)))
print(paste("Good genes accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )",sep=""))

