#functions used in MBD-seq_analysis1-5

window_plot = function(size, X, Y, dat, cut){
  windows = quantile(dat[,X], probs = seq(0, 1, by = size), na.rm = T)
  if (QUANTILE == F){
    width = max(dat[,X]) * size
    windows = seq(min(dat[,X]), max(dat[,X]), by = width)
  }
  print(windows)
  mn = c()
  x = c()
  N = c()
  sterr = c()
  cut.count = 0
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat[,X] >= left,]
    sub = sub[sub[,X] < right,]
    n = length(sub[,1])
    if (n < cut){
      cut.count = cut.count + 1
      next
    }
    mn = append(mn,mean(sub[,Y]))
    sterr = append(sterr,std.error(sub[,Y]))
    x = append(x, median(sub[,X]))
    N = append(N, n)
  }
  print(paste(cut.count, "Windows were cut Because they had below", cut, "genes"))
  plot_dat = data.frame(mn,x,sterr, N)
  return(plot_dat)
}##function to get the window data to plot th

#function to plot error bars for a set of quantiles/windows for two variables
do.window.plot = function(size, x, y, df, cut, xlab, ylab, main, loess, span, limits, xlim, ylim){
  wdat = window_plot(size, x, y, df, cut)
  if (limits == TRUE){
  	plotCI(wdat$x, wdat$mn, uiw = wdat$sterr, liw = wdat$sterr, add = F, scol = 'black', lwd = line.width, pch = 26, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, axes = F)
  	# axis(1, at = signif(seq(from = max(wdat$x), to = min(wdat$x), by = -(max(wdat$x) - min(wdat$x))/3), digits = 2))
  	# axis(1, cex.axis = cex.axis)
  	# ydif = max(wdat$mn) - min(wdat$mn)
  	# axis(2, at = signif(seq(from = min(wdat$mn), to = max(wdat$mn), by = ydif/3), digits = 2), cex.axis = cex.axis, cex.lab = cex.lab)
  } else {
  	plotCI(wdat$x, wdat$mn, uiw = wdat$sterr, liw = wdat$sterr, add = F, scol = 'black', lwd = line.width, pch = 26, main = main, xlab = xlab, ylab = ylab, axes = F)
  	#don't put in with funciton, it's easier to do case by case
  	# axis(1, at = signif(seq(from = max(wdat$x), to = min(wdat$x), by = -(max(wdat$x) - min(wdat$x))/3), digits = 2))
  	# axis(1, cex.axis = cex.axis, cex.lab = cex.lab)
  	# ydif = max(wdat$mn) - min(wdat$mn)
  	# axis(2, at = signif(seq(from = min(wdat$mn), to = max(wdat$mn), by = ydif/2), digits = 3), cex.axis = cex.axis, cex.lab = cex.lab)
  	}
  if (loess == TRUE){
    loess_fit <- loess(mn ~ x, wdat, span = span, se = T)
    lines(wdat$x, predict(loess_fit),col="red",lwd=1)
  }
  return(wdat)
}


plot.counts = function(dat, X, Y, threshold, size, xlab, ylab, span){
  Ns = c()
  Xs = c()
  windows = quantile(dat[,X], probs = seq(0, 1, by = size))
  print(windows)
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat[,X] >= left,]
    sub = sub[sub[,X] < right,]
    sig.sub = sub[sub[,Y] < threshold,]
    Ns = append(Ns, nrow(sig.sub))
    Xs = append(Xs, mean(sub[,X]))
  }
  results = data.frame(Ns, Xs)
  print(results)
  plot(Ns~Xs, data = results, xlim = range(results$Xs), xlab = xlab, ylab = ylab, axes = F)
  loess_fit <- loess(Ns~Xs, results, span = span, se = T)
  lines(results$Xs, predict(loess_fit),col="red",lwd=1)
  return(results)
}

#funciton to build randomized vectors from two columns in a dataframe
#the output can be used to repeat figures with randomized assignments for X and Y variables
sanity = function(dat, X, Y){
	x = dat[,X]
	y = dat[,Y]
	rand.x = sample(x, length(x), replace = F)
	rand.y = sample(y, length(y), replace = F)
	result = data.frame(rand.x, rand.y)
	colnames(result) = c(X, Y)
	return(result)
}


#function to perform fishers exact test for methylation classes
do.meth.fisher = function(dat, sep, cut, X, Y, alt){
	meth = na.omit(dat[dat[,X] >= sep,])
	not = na.omit(dat[dat[,X] < sep,])
	sig.meth = nrow(meth[meth[,Y] < cut,])
	other.meth = nrow(meth[meth[,Y] > cut,])
	sig.not = nrow(not[not[,Y] <= cut,])
	other.not = nrow(not[not[Y] > cut,])
	M = c(sig.meth, other.meth)
	U = c(sig.not, other.not)
	tab = as.table(rbind(M, U))
	colnames(tab) = c('sig', 'not.sig')
	rownames(tab) = c('meth', 'unmeth')
	print(tab)
	test = fisher.test(tab, alternative = alt)
	print(test)
	return(test)
}






#function to add a 'window-line' (trace of quantile means) to a plot
add.window.line = function(size, x, y, df, cut, loess, span, line.width, color){
  wdat = window_plot(size, x, y, df, cut)
  if (loess == TRUE){
    loess_fit <- loess(mn ~ x, wdat, span = span, se = T)
    lines(wdat$x, predict(loess_fit), col = color, lwd = line.width)
  }
}


#function to prepare a set of pairwise substitution rates for plotting
filter.sub.rates = function(df, sp, mut.type, meth.type, remove.zero, filter, factor, filter.mut){
  dat = df[df$species == sp,]
  if (remove.zero == TRUE){
    sub = dat[dat$dS > 0,]
  } else {
  	sub = dat
  	print(sub)
  }
  print(head(sub))
  if (filter == TRUE){	
  	#filter based on dS any gene with rate 4x higher than interquartile range
  	qnt = quantile(sub[, filter.mut], na.rm = T)
  	H = 4 * IQR(sub[, filter.mut], na.rm = T)
  	sub1 = sub[sub[, filter.mut] < qnt[2] + H,]
  } else {
  	sub1 = sub
  	print(sub)
  }
  meth.dat = sub1[,meth.type]
  mut.dat = sub1[,mut.type]
  sub2 = na.omit(data.frame(meth.dat, mut.dat))
  colnames(sub2) = c(meth.type, mut.type)
  return(sub2)
}


#function to plot just the loess lines for error bar plots
#plots multiple lines on the same panel with a selected color pallet
plot.stacked.lines = function(spp.list, color.list, xlim, ylim){
	plot.new()
	plot.window(xlim = xlim, ylim = ylim)
	axis(1, cex.axis = cex.X, line = 0, at = c(-2, 0, 2), padj = 0.25)
	axis(2, las = 1, cex.axis = 1, at = signif(seq(from = min(ylim), to = max(ylim), by = (max(ylim) - min(ylim))/3), digits = 2), cex.axis = cex.Y, line = 0)
	rho.values = c()
	p.values = c()
	for (i in 1:length(spp.list)){
		sp = spp.list[i]
		print(paste("species =", sp))
		color = color.list[i]
		sub = filter.sub.rates(dnds.dat, sp, mut.type, meth.type, remove.zeros, filter, factor, filter.mut)
		print("-----------------")
		print(head(sub))
		z = cor.test(sub[,meth.type], sub[,mut.type], method = "spearman")
		rho = z$estimate
		p = z$p.value
		rho.values = append(rho.values, rho)
		p.values = append(p.values, p)
		add.window.line(window, meth.type, mut.type, sub, 10, T, span, line.width, color)
	}
	print(paste('rhos', rho.values))
	title(main = paste(round(median(rho.values), 2), round(median(p.values), 5)))
}

#function to be used with plot.stacked.lines
#plots a set of horizontal lines labeled with species names 
#so you can double bheck which color goes with which
plot.pallet = function(spp.list, color.list){
	plot.new()
	par(mar = c(5, 15, 4, 2) + 0.1)
	plot.window(xlim = c(0,2), ylim = c(0, length(spp.list)))
	for (i in 1:length(spp.list)){
		sp = spp.list[i]
		color = color.list[i]
		abline(h = i, lwd = 8, col = color)
		
	}
	axis(2, at = 1:length(spp.list), labels = spp.list, las = 1)
	par(mar = c(5, 4, 4, 2) + 0.1)
}

#function plot.lm
#plots two variables in a dataframe and their linear regression
alpha = 0.075
library(scales)
plot.lm = function(x, y, df, xlab, ylab, main, limits, xlim, ylim){
	z = cor.test(df[,x], df[,y], method = "spearman")
	rho = signif(z$estimate, digits = 2)
	p = z$p.value
	if (limits == T){
		plot(df[,y] ~ df[,x], xlab = xlab, ylab = ylab, main = "\n", col = alpha('grey', alpha), xlim = xlim, ylim = ylim, axes = F, cex.lab = cex.lab)
		# axis(1, at = signif(seq(from = xlim[1], to = xlim[2], by = -(xlim[1] - xlim[2])/3), digits = 2))
		axis(1, cex.axis = cex.axis)
		axis(2, at = signif(seq(from = ylim[1], to = ylim[2], by = (ylim[2] - ylim[1])/2), digits = 2), cex.axis = cex.axis)
		} else {
		plot(df[,y] ~ df[,x], xlab = xlab, ylab = ylab, main = paste("rho = ", rho), col = alpha('grey', alpha))
		}
 	lm1 = lm(df[,y] ~ df[,x])
 	print(summary(lm1))
 	abline(lm1, col = 'red')
 	return(c(rho, p))
}


plot.meth.bars = function(dat, X, Y, sep, COLORS){
  m = dat[dat[,X] >= sep,]
  u = dat[dat[,X] < sep,]
  mns = c(mean(m[,Y]), mean(u[,Y]))
  meds = c(median(m[,X]), median(u[,X]))
  ses = c(std.error(m[,Y]), std.error(u[,Y]))
  print(paste("X positions:", meds))
  print(paste("Y positions:", mns))
  percent.dif = (mns[2] - mns[1]) / mns[1]
  print(paste("Percent Greater =", percent.dif))
  plotCI(meds, mns, uiw = ses, col = COLORS, add = T, pch = 19, lwd = line.width, cex = .01)
  return(percent.dif)
}

#set up lrt function that returns p values from two vectors of likelihoods
lrt = function(la, lo, df){
  G = 2*(la - lo)
  p.values = pchisq(G, df, ncp = 0, lower.tail = F)
}


#function for assign cg/gc categories of codons
get.gc.categories = function(dru){
	dru$cpg = 0
	dru$double = 0
	dru[grep('CG', dru$codons), "cpg"] <- '1'
	dru[grep('GC', dru$codons), "double"] <- '1'
	dru[grep('GG', dru$codons), "double"] <- '1'
	dru[grep('CC', dru$codons), "double"] <- '1'
	cg = na.omit(dru[dru$cpg == 1,])
	other = na.omit(dru[dru$double == 1,])
	other = other[other$cpg == 0,]
	cg = cg[order(cg$codons),]
	other = other[order(other$codons),]
	nrow(cg)
	nrow(other)
	res = rbind(cg, other)
	return(res)
}


