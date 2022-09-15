#!/bin/bash
# for nominal 1
DIR=/home/users/a/avalosma/scratch/12_TRIPLETS/not_signif/merged_nominal1
DIR2=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_all

mkdir -p $DIR2

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                file=$DIR/${data_type}_${cell_type}_CRD_gene_var_ALL.txt.gz
                echo $file
                zcat $file | awk '{ print $8";"$1 }' > $DIR2/${data_type}_${cell_type}_triplet.txt
        done
done

