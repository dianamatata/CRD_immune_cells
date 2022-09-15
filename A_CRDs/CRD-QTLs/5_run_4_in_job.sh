#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/CRD-QTLs/4_run_conditional_pass.sh"
wsbatch -J qtl.job --partition=mono-shared-EL7 --time=10:00:00 -o 4.out -e 4.err --wrap="$cmd"
