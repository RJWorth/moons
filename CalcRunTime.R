###############################################################################
### Calculate average runtime per sim
runtime=read.table('runtime.txt')

stopifnot(dim(runtime)[2]==5)

colnames(runtime)=c('Dir','machine','n','user','t')
runtime$t=runtime$t/3600
attach(runtime)

machines=levels(machine)
avg=rep(0, length(machines))
for (i in 1:length(machines))	{
		avg[i]=sum(t[machine==machines[i]])/sum(n[machine==machines[i]])	}

print('runtimes in hrs')
print(machines)
print(avg)

############# set up plotting stuff
palette(rainbow(length(machines)))
pdf('runtimes.pdf',height=8,width=4)
par(mfrow=c(2,1))
############# plot the total time vs. number of iterations
plot(n,t, pch=20, col=machine, xlab='Iterations',ylab='Time (hrs)')
for (i in 1:length(machines)) {
	abline(lm(t[machine==machines[i]]~n[machine==machines[i]]),col=i)	}
legend('topleft',pch=20,col=1:length(machines), 
	legend=machines)
############# plot time/iterations
plot(t/n,col=machine,pch=20, xlab='batch',ylab='Time per sim (hrs)')
for (i in 1:length(machines)) abline(h=avg[i],col=i)
############# finish plot
dev.off()

detach(runtime)


