#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/6.4_run_mapping_CRD_gene_conditional.sh"
wsbatch -J 6.4script.job --partition=mono-EL7 --time=04:00:00 -o 6.4script.out -e 6.4script.err --wrap="$cmd"


