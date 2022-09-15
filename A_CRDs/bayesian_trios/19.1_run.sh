#!/bin/bash

for data_type in  'methyl' 'hist' ; do
        for cell1 in 'neut' 'mono' 'tcell' ; do
		for cell2 in 'neut' 'mono' 'tcell' ; do
			if [ "$cell1" != "$cell2" ]; then
				echo $cell1 $cell2 $data_type
				file=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/bayesian_trios/19_triplets_sharing_btw_cells.R
				cmd="Rscript $file $data_type $cell1 $cell2"
				wsbatch -J 19.job --partition=shared-bigmem --time=12:00:00 --mem=500000 -o 19.out -e 19.err --wrap="$cmd"
			fi
        	done
	done
done
