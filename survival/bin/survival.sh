#!/bin/bash
#./survival single/pairs file outfile
mode=$1
file=$2
outFile=$3

if [ $mode == 'single' ]
then
	#for i in `awk '{if($7<=0.05) print $1}' GenesSinglesScoreFisherSorted.txt`;
	while read -r a ; do ./callLifelines.sh $a $3; done <$2
	#do
		#./callLifelines.sh $i;
	#done
elif [ $mode == 'pairs' ]
then
	#awk '{if($10<=0.05 && $12<=0.05) {split($1,a,"|"); print a[1],a[2]; }}' GenesPairsScoreFisher.txt > pairGenes_significant_both
	while read -r a b ; do ./callLifelines.sh $a $b $3; done <$2
fi
