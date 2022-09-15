#!/bin/bash

DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/A_CRD_plots

for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
        cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/1c_plot_mixture_model_methylCRDs.R '$cell_type' $DATADIR $OUT_FOLDER"
	eval $cmd
done

#'methyl_neut' /home/users/a/avalosma/scratch/1_CRD /home/users/a/avalosma/scratch/A_CRD_plots
