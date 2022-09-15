#!/bin/bash

DATADIR=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
OUTDIR=$DATADIR/mapping_gene_CRDs
mkdir $OUTDIR/merged

for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                        name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
                        for module in 'mean' 'loom' ; do
#                        	for condition in 'mapping_CRD_gene' 'inverse_mapping_CRD_gene' 'mapping_CRD_gene_nominal' 'inverse_mapping_CRD_gene_nominal'; do
				echo $name
				condition='mapping_CRD_gene'
                                cat $OUTDIR/${module}/${condition}/${name}_*.txt | gzip -c > $OUTDIR/merged/${name}_${module}_${condition}_permuts.txt.gz
#                       	done
			done
                done
        done
done

