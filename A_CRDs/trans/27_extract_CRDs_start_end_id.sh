#!/bin/bash


DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/7_CRD_Trans/CRD_id_start_end
mkdir -p $OUTDIR

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		module='mean'
		MOD=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
		zcat $MOD | cut -f2-4 > $OUTDIR/${data_type}_${cell_type}_CRD_id.txt
	done
done

