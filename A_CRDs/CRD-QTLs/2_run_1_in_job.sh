#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/CRD-QTLs/1_run_permutation_pass.sh"
wsbatch -J qtl.job --partition=bigmem-EL7 --time=36:50:00 -o OUT/1p.out -e OUT/1p.err --wrap="$cmd"

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/CRD-QTLs/1_run_nominal_pass.sh"
wsbatch -J qtlnom.job --partition=bigmem-EL7 --time=36:50:00 -o OUT/1n.out -e OUT/1n.err --wrap="$cmd"

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/CRD-QTLs/1_run_permutation_pass.sh"
wsbatch -J qtl1000.job --partition=shared-bigmem-EL7 --time=12:00:00 -o OUT/1p1000.out -e OUT/1p1000.err --wrap="$cmd"



cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/CRD-QTLs/1_run_nominal_pass_no_FDR.sh"
wsbatch -J qtl_nom_noFDR.job --partition=shared-bigmem-EL7 --time=12:00:00 -o OUT/qtl_nom_noFDR.job.out -e OUT/qtl_nom_noFDR.job.err --wrap="$cmd"