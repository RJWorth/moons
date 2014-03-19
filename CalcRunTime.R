###############################################################################
### Calculate average runtime per sim
runtime=read.table('runtime.txt')
#runtime=read.table(paste(prefix,'runtime.txt',sep=''))

if(dim(runtime)[2]==5) colnames(runtime)=c('Dir','machine','n','user','t')
#if(dim(runtime)[2]==2) {
#	colnames(runtime)=c('Dir','t')
#	machine=rep('nova',length(runtime$Dir))
#		machine[grep('MDir',runtime$Dir)]='myra'
#		machine[grep('CDir',runtime$Dir)]='chloe'
#		machine=as.factor(machine)
#		n=rep(1,length(runtime$Dir))
#		runtime=cbind(runtime,machine,n)	}
runtime$t=runtime$t/3600
attach(runtime)

machines=levels(machine)
avg=rep(0, length(machines))
for (i in 1:length(machines))	{
		avg[i]=sum(t[machine==machines[i]])/sum(n[machine==machines[i]])	}

print('runtimes in hrs')
print(machines)
print(avg)

palette(rainbow(length(machines)))
pdf('runtimes.pdf',height=8,width=4)
par(mfrow=c(2,1))
##########
plot(n,t, pch=20, col=machine, xlab='Iterations',ylab='Time (hrs)')
for (i in 1:length(machines)) {
	abline(lm(t[machine==machines[i]]~n[machine==machines[i]]),col=i)	}
legend('topleft',pch=20,col=1:length(machines), 
	legend=machines)
#############
plot(t/n,col=machine,pch=20, xlab='batch',ylab='Time per sim (hrs)')
for (i in 1:length(machines)) abline(h=avg[i],col=i)

dev.off()

detach(runtime)

