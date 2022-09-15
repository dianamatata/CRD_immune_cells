#!/bin/bash
DIR_T=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_data
DIR_BN=/home/users/a/avalosma/scratch/12_TRIPLETS/BN
mkdir -p $DIR_BN

for file in $DIR_T/*.txt ; do
	echo $file
	f=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'_' -f1-3)
	Rscript 16_bnlearn.call2.R $file $DIR_BN/BN_$f
done
