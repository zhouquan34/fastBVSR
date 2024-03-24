
args<-commandArgs(TRUE)
dat = read.table(paste(args[1], '.path.txt', sep=''), header=T)
png(filename=paste(args[1],'_trace.png', sep=''),res=72,unit='in',width=12,height=18)
par(mfrow=c(3,1))
par(mar=c(6,6,2,2))
plot(dat[,c(1,3)],xlab='Iteration', ylab='Model size', type='l', cex.lab=1.5,cex.axis=1.5)
plot(dat[,c(1,4)],xlab='Iteration', ylab='Heritability', type='l',cex.lab=1.5,cex.axis=1.5)
plot(dat[,c(1,5)],xlab='Iteration', ylab='PVE', type='l',cex.lab=1.5,cex.axis=1.5)
dev.off()

