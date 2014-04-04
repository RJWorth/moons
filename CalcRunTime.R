###############################################################################
### Calculate average runtime per sim
runtime=read.table('runtime.txt')

stopifnot(dim(runtime)[2]==6)

colnames(runtime)=c('Dir','machine','iter','n','user','t')
attach(runtime)

ntot=n*iter	#number of rocks integrated, including via iteration
machines=levels(machine)
avg=rep(0, length(machines))
for (i in 1:length(machines))	{
		avg[i]=sum(t[machine==machines[i]])/sum(iter[machine==machines[i]]*n[machine==machines[i]])	}

print('runtimes in s')
print(machines)
print(avg)

############# set up plotting stuff
palette(rainbow(length(machines)))
pdf('runtimes.pdf',height=4,width=8)
par(mfrow=c(1,2))
############# plot 1
plot(iter,t/(255+n), pch=20, col=machine)
	fit2=matrix(data=NA,nrow=length(machines),ncol=2)
for (i in 1:length(machines)) {
	fit2[i,]=lm(t[machine==machines[i]]/(255+n[machine==machines[i]])~iter[machine==machines[i]])$coefficients
	abline(fit2[i,1],fit2[i,2],col=i)	}
legend('bottomright',pch=20,col=1:length(machines), 
	legend=machines)
############# plot 2
plot(n, t/iter, col=machine,pch=20)
	fit1=matrix(data=NA,nrow=length(machines),ncol=2)
for (i in 1:length(machines)) {
#for (i in c(1,3,4,5)) {
	fit1[i,]=lm(t[machine==machines[i]]/iter[machine==machines[i]]~n[machine==machines[i]])$coefficients
	abline(fit1[i,1],fit1[i,2],col=i)	}
colnames(fit1)=c('LatencyPerSim','TimePerObject')
rownames(fit1)=machines
print(fit1)
############# finish plot
dev.off()

detach(runtime)
############### function to predict runtime
RunTime=function(n, iter, mach, fit1)	{
	latency=fit1[machines==mach, 1]
	slope=fit1[machines==mach, 2]

	t=(latency + slope*n)*iter
	print(paste("Time in hrs =",t/3600))
	return(t)	}

############### function to give max n for given timeframe (hrs)
NObjects=function(t, iter, mach, fit1)	{
	latency=fit1[machines==mach, 1]
	slope=fit1[machines==mach, 2]

	n=floor(((t*3600)/iter-latency)/slope)
	print(paste('Max n for',iter,'iterations,',t,'hrs =',n))
	return(n)	}

############### 




