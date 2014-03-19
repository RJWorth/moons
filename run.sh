#!/bin/sh
############################################################################### 
### Script to repeatedly run moon-planet collision ratio simulations
# Start time
t1=$(date +%s)
machine=$(hostname -s)

### Simulation parameters
time=3		# = log(years)
output=1	# = log(years)
step=0.5	# = days
niter=10	# = number of iterations to run

vers='mercury_TidesGas.for'
user='no'	# use user-defined forces?

nobj=1000	# number of ejected fragments
pl='J'		# planet to aim for
mode='gen'	# 'gen' to generate rock list, 'good' to use it

### clean?
#\rm $1/infosum.out
#\rm runtime.txt

### Write files.in
./writefiles.sh $1

### Repeat simulation n times and record summary of collisions
### Range for iterations
if [ $machine = chloe ]; then
	itrange=$(jot $niter 1)	# start at y, go up x-1 steps
else
	itrange=$(seq 1 $niter)	# go from x to y
fi

### Do iterations
	for j in $itrange
	do
	echo 'run: '$1', '$j
	\rm $1/Out/*.dmp
	\rm $1/Out/*.out

#### Randomize moon phases and rock positions
		if [ $mode = 'gen' ]; then
		python -c 'import MercModule; MercModule.makebigrand("'$1'","'$j'")'
		python -c 'import MercModule; MercModule.makesmall("'$1'","'$j'","'$nobj'","'$pl'",1.0e12,1.0e2)'
#### OR choose moon phases, use rocks from good.in
		elif [ $mode = 'good' ]; then
		python -c 'import MercModule; MercModule.makebigchoose("'$1'","'$j'")'
		python -c 'import MercModule; MercModule.good2small("'$1'","'$j'","200")'
		fi
	# Write param.in file
	./writeparam.sh $1 $time $output $step $time $user
	# Compile mercury
	gfortran -w -o $1/Out/merc_$1 Files/$vers

#### Run mercury
		cd $1/Out;	./merc_$1;	cd ../..

### Write collisions summary, copy good in coords
		if [ $mode = 'gen' ]; then
		python -c 'import MercModule; MercModule.copyinfo("'$1'","'$j'",True)'
		elif [ $mode = 'good' ]; then
		python -c 'import MercModule; MercModule.copyinfo("'$1'","'$j'",False)'
		fi
	done	# j iterations

### Send alert to my email
	./email.sh $1 $niter 'Moons'
# Write stop time for this directory:
t2=$(date +%s)
echo $1"	"$machine"	"$niter"	"$user"	"$(echo "$t2 - $t1"|bc ) >> runtime.txt

