#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/1_quantify_CRD.sh"
wsbatch -J B_2_script.job --partition=mono-EL7 --time=01:10:00 -o B_2_script.out -e B_2_script.err --wrap="$cmd"
