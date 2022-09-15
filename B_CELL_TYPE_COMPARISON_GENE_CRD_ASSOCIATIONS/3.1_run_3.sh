#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/3_run_mapping_CRD_gene.sh"
name="B3_script"
wsbatch -J ${name}.job --partition=mono-EL7 --time=10:10:00 -o ${name}.out -e ${name}.err --wrap="$cmd"
