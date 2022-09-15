#!/bin/bash

# CLOMICS build to make the dendogram/tree, topo to annotate and call to have a threshold
# next filter on col 25 N_REG col 30 MOD, MOD=1 if CRD, N_REG nbr of peaks, at least 2 RE

# Olivier's code for CRD
CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics

# Input data folder so far because I haven't generated yet
DATADIR=/home/users/a/avalosma/scratch/0_CRD

# My output folder for CRD 2nd part of code project
OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
        LI=$DATADIR/${cell_type}.bed.gz
        for c in $(seq 1 22); do
                LO=$OUT_FOLDER/$cell_type\.chr$c
                echo "$cell_type $c"
                cmd1="$CLOMICs build --bed $LI --region $c --out $LO.1.tgz --silent"
                cmd2="$CLOMICs topo --bed $LI --tree $LO.1.tgz --chr $c --out ${LO}.tree.txt.gz && rm ${LO}.*.tgz"
                cmd3="$CLOMICs call --tree ${LO}.tree.txt.gz --threshold 2 --out ${LO}.module.txt.gz"
                cmd4="bgzip ${LO}.module.txt.gz"
                cmd="$cmd1 && $cmd2 && $cmd3"
#                echo "$cmd"
                FILENAME=$OUT_FOLDER/OUT/${cell_type}_${c}_clomics
        	wsbatch -J $FILENAME.job --partition=shared-EL7 --time=10:00:00 -o $FILENAME.out -e $FILENAME.err --wrap="$cmd"
	done
done

for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
        LDIST=$OUT_FOLDER/${cell_type}_dist.ALLchr.bed
        for c in $(seq 1 22); do
                LT=$OUT_FOLDER/${cell_type}.chr$c\.module.txt.gz
                zcat $LT | awk -v OFS='\t' -v chr=$c '{if($30 == 1) print $19}' >> $LDIST
        done
done
