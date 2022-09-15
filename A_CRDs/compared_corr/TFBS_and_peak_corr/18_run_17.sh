#!/bin/bash
FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
for pair in "mono_neut" "mono_tcell" "neut_tcell"
    do
	cmd="Rscript $FOLDER/TFBS_and_peak_corr/17_plot_heatmap_per_pair.R"
	echo $cmd
	wsbatch -J odd_${pair} \
	--partition=shared-cpu \
	--time=12:00:00 \
	--mem=10000 \
	-o $FOLDER/log/odd_${pair}.out \
	-e $FOLDER/log/odd_${pair}.err \
	--wrap="$cmd"
done
