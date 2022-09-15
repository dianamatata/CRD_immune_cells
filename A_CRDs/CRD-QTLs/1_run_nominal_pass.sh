#!/bin/bash

VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs
mkdir -p $OUTDIR $OUTDIR/nominal $OUTDIR/merged_nominal
K=100
FDR=0.05

for data_type in  'methyl' 'hist' ; do
  for cell_type in 'neut' 'mono' 'tcell' ; do
    for module in 'mean' 'loom' ; do
      BED=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
      for k in $(seq 1 $K); do
        OUT1=$OUTDIR/nominal/${data_type}_${cell_type}_${module}_CRD_QTL_$k
        cmd="QTLtools cis --vcf $VCF --bed $BED --nominal $FDR  --chunk $k $K --out ${OUT1}.txt"
        eval $cmd
      done
    done
  done
done

### concatenate at the end

for data_type in  'methyl' 'hist' ; do
  for cell_type in 'neut' 'mono' 'tcell' ; do
    for module in 'mean' 'loom' ; do
      name=${data_type}_${cell_type}_${module}
      cat $OUTDIR/nominal/${name}_*.txt | gzip -c > $OUTDIR/merged_nominal/${name}_${condition}_ALL.txt.gz
    done
  done
done
