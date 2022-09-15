#!/bin/bash

DATADIR=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
mkdir -p $DATADIR/quantify_ALL

for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                        name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
                        for module in 'mean' 'loom' ; do
 
                                # unzip all the mean.txt.gz files
                                gunzip $DATADIR/quantify/${name}.chr*.${module}.txt.gz
 
                                # take the first file fully and then everything except the first line into unsorted
                                awk 'FNR>1||NR==1' $DATADIR/quantify/${name}.chr*.${module}.txt \
				 > $DATADIR/quantify_ALL/${name}.ALLchr.${module}.unsorted.txt
 
                                #sort all the lines except header, by col 1 first and then col 2 by n numerical value
                                awk 'NR == 1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' \
				$DATADIR/quantify_ALL/${name}.ALLchr.${module}.unsorted.txt | bgzip -c \
				> $DATADIR/quantify_ALL/${name}.ALLchr.${module}.txt.gz
 
                                tabix -p bed $DATADIR/quantify_ALL/${name}.ALLchr.${module}.txt.gz
                                bgzip $DATADIR/quantify/${name}.chr*.${module}.txt
                                rm  $DATADIR/quantify_ALL/${name}.ALLchr.${module}.unsorted.txt
                        done
                done
        done
done
