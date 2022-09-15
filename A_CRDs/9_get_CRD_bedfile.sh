#!/bin/bash

# This script takes the CRD ID (ex: H3K4me1_chr10-51569), the start and end point of the CRD, and stack all the results of all the chromosomes in 1 file

DATADIR2=/home/users/a/avalosma/scratch/1_CRD
DATADIR3=/home/users/a/avalosma/scratch/3_CRD

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	for c in $(seq 1 22); do
        	LT=$DATADIR2/${cell_type}.chr$c\.module.txt.gz
		LOUT=$DATADIR3/${cell_type}.chr$c\.module.bed
		zcat $LT | awk -v OFS='\t' -v chr=$c '{if($30==1 && $25>1) print chr,$17,$18,$4}' > $LOUT
		cat $DATADIR3/${cell_type}.chr*.module.bed | sort -V -k1,1 -k2,2n > $DATADIR3/${cell_type}.ALLchr.module.bed
	done
done


