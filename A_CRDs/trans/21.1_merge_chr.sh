#!/bin/bash
OUTDIR=/home/users/a/avalosma/scratch/7_CRD_Trans

### concatenate at the end
mkdir -p $OUTDIR/merged

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
			OUT=$OUTDIR/OUTsimple/${data_type}_${cell_type}_${module}_chr_${c0}_${c1}_trans
			cat $OUTDIR/OUTsimple/${data_type}_${cell_type}_${module}_* \
			> $OUTDIR/merged/${data_type}_${cell_type}_${module}_trans_ALL.txt.gz
                done
        done
done
