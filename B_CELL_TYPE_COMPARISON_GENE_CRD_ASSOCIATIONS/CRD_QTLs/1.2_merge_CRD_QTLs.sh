#!/bin/bash
permuts=1000
K=100

OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs


for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                        name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
                        for module in 'mean' 'loom' ; do
                                cat $OUTDIR/permuts_$permuts/${name}_${module}_CRD_QTL_*.txt | gzip -c > $OUTDIR/merged_$permuts/${name}_${module}_permuts.txt.gz
			done
		done
	done
done
