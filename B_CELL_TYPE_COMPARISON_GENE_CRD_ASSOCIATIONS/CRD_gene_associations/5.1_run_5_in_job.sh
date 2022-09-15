#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/5_mapping_CRD_gene_conditional.sh"
wsbatch -J condi.job --partition=mono-EL7 --time=04:00:00  -e condi.err --wrap="$cmd"
