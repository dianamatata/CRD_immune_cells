#!/bin/bash
permuts=1000
K=100

VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
DATADIR=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs
mkdir -p $OUTDIR $OUTDIR/permuts_$permuts $OUTDIR/merged_$permuts

for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                        name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
			for module in 'mean' 'loom' ; do
                                MOD=$DATADIR/quantify_ALL/${name}.ALLchr.${module}.txt.gz
				for k in $(seq 1 $K); do
					OUT1=$OUTDIR/permuts_$permuts/${name}_${module}_CRD_QTL_$k
					cmd="QTLtools cis --vcf $VCF --bed $MOD --permute $permuts  --chunk $k $K --out ${OUT1}.txt"
					eval $cmd
				done
			done
		done
	done
done
