#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/5_run_mapping_CRD_gene_permutation_pass.sh"
wsbatch -J 5script.job --partition=mono-EL7 --time=04:00:00 -o 5script.out -e 5script.err --wrap="$cmd"


