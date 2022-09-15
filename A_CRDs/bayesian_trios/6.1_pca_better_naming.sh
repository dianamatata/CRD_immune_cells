#!/bin/bash

DIR_2col=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1000/2_columns
DIR_PC=/home/users/a/avalosma/scratch/12_TRIPLETS/PC1

# for nominal 1 ie ALL
DIR_2col=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1/2_columns
DIR_PC=/home/users/a/avalosma/scratch/12_TRIPLETS/PC1_all

## merge 2 columns in 1

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		awk '{ print $1";"$2 }' $DIR_2col/${data_type}_${cell_type}_mean_mapping_CRD_gene_ALL.txt \
		> $DIR_2col/${data_type}_${cell_type}_mean_mapping_CRD_gene_merged.txt
	done
done


## replace this new column with col 4
# file 2 single column, replace col in file 1
col=4
for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		file1=$DIR_PC/quantify/${data_type}_${cell_type}.pc1.txt
		file2=$DIR_2col/${data_type}_${cell_type}_mean_mapping_CRD_gene_merged.txt
		awk 'FNR==NR{a[NR]=$0;next} {sub($4, a[FNR])}1' $file2 $file1 \
		> $DIR_PC/quantify/${data_type}_${cell_type}.pc1.new.txt
	done
done

## bgzip and tabix

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                OUTPUT=$DIR_PC/quantify/${data_type}_${cell_type}
                awk 'NR == 1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' $OUTPUT.pc1.new.txt \
                 | bgzip -c > $OUTPUT.pc1.sorted.txt.gz
                tabix -p bed $OUTPUT.pc1.sorted.txt.gz
        done
done
