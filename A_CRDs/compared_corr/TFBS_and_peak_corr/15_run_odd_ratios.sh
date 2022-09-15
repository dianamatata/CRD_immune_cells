#!/bin/bash
FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
for chr in {1..22}
do
        for pair in "mono_neut" "mono_tcell" "neut_tcell"
        do
		cmd="Rscript $FOLDER/TFBS_and_peak_corr/14_compute_odd_ratios_per_chr_pair_sep.R $chr $pair"
		echo $cmd
		wsbatch -J odd_${chr}_${pair} \
		--partition=shared-cpu \
		--time=12:00:00 \
		--mem=10000 \
		-o $FOLDER/log/odd_${chr}_${pair}.out \
		-e $FOLDER/log/odd_${chr}_${pair}.err \
		--wrap="$cmd"
	done
done
