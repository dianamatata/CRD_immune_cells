#!/bin/bash

DATADIR=/home/users/a/avalosma/scratch/2_CRD

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
        	for module in 'mean' 'loom' ; do
			# unzip all the mean.txt.gz files
			gunzip $DATADIR/quantify/${data_type}_${cell_type}.chr*.${module}.txt.gz
	 
			# take the first file fully and then everything except the first line into unsorted
			awk 'FNR>1||NR==1' $DATADIR/quantify/${data_type}_${cell_type}.chr*.${module}.txt \
			  > $DATADIR/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.unsorted.txt
			
			#sort all the lines except header, by col 1 first and then col 2 by n numerical value
			awk 'NR == 1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' \
			 $DATADIR/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.unsorted.txt \
			 | bgzip -c > $DATADIR/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
	 
			tabix -p bed $DATADIR/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
			bgzip $DATADIR/quantify/${data_type}_${cell_type}.chr*.${module}.txt
                	rm  $DATADIR/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.unsorted.txt
		done        
	done
done

for file in $DATADIR/quantify/*; do gzip "$file"; done

