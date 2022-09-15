#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/CRD_QTLs/1bis_run_permutation_pass.sh"

# time: days-hours:minutes:seconds
wsbatch -J qtl1000.job --partition=public-cpu --time=4-00:00:00 -o 1p1000pu.out -e 1p1000pu.err --wrap="$cmd"
# time limit problem with:
# wsbatch -J qtl1000.job --partition=shared-bigmem --time=12:00:00 -o 1p1000.out -e 1p1000.err --wrap="$cmd"
# (ReqNodeNotAvail, Reserved for maintenance)
# wsbatch -J qtl1000.job --partition=public-bigmem --time=4-00:00:00 -o 1p1000.out -e 1p1000.err --wrap="$cmd"


cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/CRD_QTLs/6.2_variants_linked_to_CRD_in_LD.sh"

wsbatch -J CRD_QTL_LD.job --partition=public-cpu --time=4-00:00:00 -o OUT/CRD_QTL_LD.out -e OUT/CRD_QTL_LD.err --wrap="$cmd"
