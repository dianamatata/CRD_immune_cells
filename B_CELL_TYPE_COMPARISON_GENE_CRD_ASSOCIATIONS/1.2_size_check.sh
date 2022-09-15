#!/bin/bash

OUT_FOLDER=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify
ANALYSIS_FILE=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/1.2b_analysis_sample_size.txt


echo '' > $ANALYSIS_FILE
for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                        name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
                        for mod in 'mean' 'loom' ; do
				count=$(cd $OUT_FOLDER | ls | grep $name | grep $mod | wc -l)
				echo $name $mod $count >>  $ANALYSIS_FILE
			done
		done
	done
done


