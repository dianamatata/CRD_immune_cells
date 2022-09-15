#!/bin/bash
# can be added as cmd 2 in file 8? 
DATADIR2=/home/users/a/avalosma/scratch/1_CRD
DATADIR3=/home/users/a/avalosma/scratch/3_CRD

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	for c in $(seq 1 22); do
        	LT=$DATADIR2/${cell_type}_94s.chr$c\.module.txt.gz
		LOUT=$DATADIR3/${cell_type}_94s.chr$c\.module.bed
		zcat $LT | awk -v OFS='\t' -v chr=$c '{if($30==1 && $25>1) print chr,$17,$18,$4}' > $LOUT
		cat $DATADIR3/${cell_type}_94s.chr*.module.bed | sort -V -k1,1 -k2,2n > $DATADIR3/${cell_type}_94s.ALLchr.module.bed
	done
done


