#!/bin/bash

                                
cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots"
                                
wsbatch -J plotCRDs.job --partition=shared-bigmem --time=12:00:00 --mem=500000 -o plotCRDs.out -e plotCRDs.err --wrap="$cmd"
