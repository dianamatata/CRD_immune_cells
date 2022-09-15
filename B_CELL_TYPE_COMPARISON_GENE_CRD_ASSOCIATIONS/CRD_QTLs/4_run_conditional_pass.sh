#!/bin/bash
permuts=1000
K=100

VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
DATADIR=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs
mkdir -p $OUTDIR $OUTDIR/conditional

for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                        name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
			for module in 'mean' 'loom' ; do
	                        TH=$OUTDIR/significants/FDR_0.05_${name}_${module}_permuts.thresholds.txt
                                MOD=$DATADIR/quantify_ALL/${name}.ALLchr.${module}.txt.gz
				for k in $(seq 1 $K); do
					OUT1=$OUTDIR/conditional/${name}_${module}_CRD_QTL_$k
					cmd="QTLtools cis --vcf $VCF --bed $MOD --mapping $TH  --chunk $k $K --out ${OUT1}.txt"
					wsbatch -J qtl1000.job --partition=shared-bigmem --time=1:00:00 --wrap="$cmd"
				done
			done
		done
	done
done

### concatenate at the end
mkdir -p $OUTDIR/conditional_merged

for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                	for module in 'mean' 'loom' ; do
                        	name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}_${module}
                        	cat $OUTDIR/conditional/${name}_*.txt | gzip -c > $OUTDIR/conditional_merged/${name}.txt.gz
			done
                done
        done
done

