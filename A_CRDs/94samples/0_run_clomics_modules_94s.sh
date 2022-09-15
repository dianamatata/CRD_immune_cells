#!/bin/bash
 
# Olivier's code for CRD
CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
 
# Input data folder so far because I haven't generated yet
DATADIR=/home/users/a/avalosma/scratch/0_CRD/94_samples
 
# My output folder for CRD 2nd part of code project
OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT
 
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	LI=$DATADIR/${cell_type}_merged_residuals_94s.bed.gz # only input	
	for c in $(seq 1 22); do
		LO=$OUT_FOLDER/${cell_type}_94s.chr$c
                echo "$cell_type $c"
		cmd1="$CLOMICs build --bed $LI --region $c --out $LO.1.tgz --silent"
		cmd2="$CLOMICs topo --bed $LI --tree $LO.1.tgz --chr $c --out $LO\.tree.txt.gz && rm $LO\.*.tgz"
		cmd3="$CLOMICs call --tree $LO\.tree.txt.gz --threshold 2 --out $LO\.module.txt.gz"
		cmd4="bgzip $LO\.module.txt.gz"
		cmd="$cmd1 && $cmd2 && $cmd3"

		# CLOMICS build to make the dendogram/tree, topo to annotate and call to have a threshold
		# next filter on col 25 N_REG col 30 MOD, MOD=1 if CRD, N_REG nbr of peaks, at least 2 RE

		# command in cluster
		echo "$cmd"
		FILENAME=$OUT_FOLDER/OUT/${cell_type}_${c}_clomics
		sbatch -J $FILENAME.job --partition=mono-EL7 --time=00:08:00 -o $FILENAME.out -e $FILENAME.err --wrap="$cmd"
	done

done
