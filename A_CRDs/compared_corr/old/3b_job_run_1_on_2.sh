#! /bin/bash

FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
cmd="source $FOLDER/old/3_run_1_on_2.sh"
wsbatch -J v3.job --partition=shared-bigmem --time=12:00:00 --mem=100000 -o $FOLDER/log/v3.out -e $FOLDER/log/v3.err --wrap="$cmd"


file="/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/files/mono_neut_tcell.chr"
FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
for chr in {1..22}; do 
	cmd="source $FOLDER/old/3_run_1_on_2_in_filename.sh $file${chr}.txt"
	echo $cmd
	wsbatch -J j_$chr --partition=shared-bigmem --time=12:00:00 --mem=100000 -o $FOLDER/log/j_$chr.out -e $FOLDER/log/j_$chr.err --wrap="$cmd"
done
