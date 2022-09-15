#!/bin/bash

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/10_correlation_peaks_vs_PCHiC.R "
#wsbatch -J 10_1.job --partition=shared-bigmem -o 10_1.out -e 10_1.err --wrap="$cmd"
wsbatch -J 10_1.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o 10_1.out -e 10_1.err --wrap="$cmd"

