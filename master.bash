#!/bin/bash
############################################################################### 
# Start an instance of runmoons.sh in each Dir

machine=$(hostname -s)
if [ $machine = basil ]; then
	Dirs=(BDir1)
elif [ $machine = chloe ]; then
	Dirs=(DirTest1)
elif [ $machine = drake ]; then
	Dirs=(DDir1)
elif [ $machine = halley ]; then
	Dirs=(HDir1)
elif [ $machine = lemay ]; then
	Dirs=(LDir1)
fi

for i in ${Dirs[*]}
do
	
	nice -n 10 ./runmoons.sh $i >& $i/run.pipe &
	echo 'master: '$i'  '$!
done

#  after above are finished, run finish.sh for summary.out file and R plots
