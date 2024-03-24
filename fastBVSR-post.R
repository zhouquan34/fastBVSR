
args<-commandArgs(TRUE)
w = as.numeric(args[1])
dat = numeric(0)
last = 0
for (i in 2:length(args)){
	d = read.table(args[i], header=T)
	k = which(d[,1] > w)
	if (length(k) > 0) {
		l = d[k[1], 1] - w 
		if (l < d[k[1],2]){
			d[k[1],2] = l
		}
		last = d[nrow(d), 1]
		d = d[k,2:5]
		dat = rbind(dat, d)
	}
}

n = last - w
qs = c(0.025, 0.05, 0.5, 0.95, 0.975)
stat = c("ModelSize", "Heritability", "PVE")
for (i in 1:3){
	qi = 1
	ave = sum(dat[,i+1] * dat[,1])/n
	ord = order(dat[,i+1])
	s = 0
	quant = numeric(5)
	for (j in 1:length(ord)){
		s = s + dat[ord[j],1]
		while (s / n >= qs[qi]){
		#	cat(i, qi, s, dat[ord[j], i+1], "\n")
			quant[qi] = dat[ord[j], i+1]
			qi = qi + 1
			if (qi == 6){break}
		}
		if (qi == 6){
			break	
		}
	}
	cat(stat[i], ":\n", sep="")
	cat("Mean = ", ave, "; Median = ", quant[3], "\n", sep="")
	cat("90% credible interval = (", quant[2], ", " , quant[4], ")\n", sep="")
	cat("95% credible interval = (", quant[1], ", " , quant[5], ")\n", sep="")
	cat("\n")
}



