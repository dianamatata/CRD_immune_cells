#!/bin/bash

DATADIR3=/home/users/a/avalosma/scratch/3_CRD

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	zcat $DATADIR3/${cell_type}.chr*.peak.txt.gz | sort -V -k1,1 > $DATADIR3/${cell_type}.ALLchr.peaksID
done

# cut  -d' ' -f1
