#!/bin/bash

DIR_BN=/home/users/a/avalosma/scratch/12_TRIPLETS/BN

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		cat BN_${data_type}_${cell_type}_* >> BN_${data_type}_${cell_type}_ALL.txt
		echo $(ls *${data_type}_${cell_type}* | wc -l)
	done
done

