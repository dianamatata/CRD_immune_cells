#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/trans/21_trans_analysis.sh"
job="21trans"
wsbatch -J ${job}.job --partition=public-cpu --time=24:00:00 -o ${job}.out -e ${job}.err --wrap="$cmd"
