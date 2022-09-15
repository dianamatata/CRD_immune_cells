#!/bin/bash

# filter the CRDs
# col 25 N_REG col 30 MOD, MOD=1 if CRD, N_REG nbr of peaks, at least 2 R

DATADIR=/home/users/a/avalosma/scratch/0_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD

for cell_type in 'neut'  'mono' 'tcell' ; do        
	LI=$DATADIR/quantif_M_hist_${cell_type}.bed.gz
        for c in $(seq 1 22); do
		LO=$OUT_FOLDER/hist_${cell_type}.chr$c
		zcat ${LO}.module.txt.gz | awk '{ if ($30 == 1 && $25 > 1) print $4}'
	done | gzip -c > $OUT_FOLDER/hist_$cell_type\_ALL.CRDs.txt.gz
done

