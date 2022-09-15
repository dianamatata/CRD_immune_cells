#!/bin/bash

OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD

for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
        LDIST=$OUT_FOLDER/${cell_type}_dist.ALLchr.bed
        for c in $(seq 1 22); do
                LT=$OUT_FOLDER/${cell_type}.chr$c\.module.txt.gz
                zcat $LT | awk -v OFS='\t' -v chr=$c '{if($30 == 1) print $19}' >> $LDIST
        done
done
