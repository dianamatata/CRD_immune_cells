#!/bin/bash

VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs
mkdir -p $OUTDIR $OUTDIR/conditional
K=100

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
                        BED=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
			TH=$OUTDIR/significants/FDR_0.05_${data_type}_${cell_type}_${module}__ALL.thresholds.txt
                        for k in $(seq 1 $K); do
                                OUT1=$OUTDIR/conditional/${data_type}_${cell_type}_${module}_CRD_QTL_$k
                                cmd="QTLtools cis --vcf $VCF --bed $BED --mapping $TH --chunk $k $K --out ${OUT1}.txt"
                                eval $cmd
                        done
                done
        done
done

### concatenate at the end
mkdir -p $OUTDIR/conditional_merged

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
                        name=${data_type}_${cell_type}_${module}
                        cat $OUTDIR/conditional/${name}_*.txt | gzip -c > $OUTDIR/conditional_merged/${name}_${condition}_ALL.txt.gz
                done
        done
done
