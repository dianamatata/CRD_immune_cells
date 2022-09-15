#!/bin/bash
OUTDIR=/home/users/a/avalosma/scratch/5_CRDgene
analysis_file=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/5.3_analysis_chunks.txt
echo ' ' > $analysis_file

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
			name=${data_type}_${cell_type}_${module}
			cd $OUTDIR/mapping_CRD_gene_permut200
			nbr_chunks=$(ls | grep $name | wc -l)
			echo $name $condition $nbr_chunks >> $analysis_file	
                done
        done
done
