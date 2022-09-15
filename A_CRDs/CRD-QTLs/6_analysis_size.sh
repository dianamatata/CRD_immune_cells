#!/bin/bash

DIR1=/home/users/a/avalosma/scratch/10_CRD_QTLs/significants
DIR2=/home/users/a/avalosma/scratch/10_CRD_QTLs/conditional_merged

# DIR1 after 5% FDR
# DIR2 after conditional pass

analysis_file=7_signif_analysis.txt
analysis_file1=7_signif_analysis_FDR5.txt

echo ' ' > $analysis_file1
for file in $DIR1/*significant.txt; do
        f=$(echo $file | rev | cut -d "/" -f1 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file1
done

echo ' ' > $analysis_file
for file in $DIR2/*.txt*; do
        f=$(echo $file | rev | cut -d "/" -f1-2 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file
done


# compare with G's values
DIRG=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/

analysis_file2=7.1_signif_analysis_Guillaume.txt

echo ' ' > $analysis_file2

for file in $DIRG/*/permute/*/*_permutations* ; do
        f=$(echo $file | rev | cut -d "/" -f1-2 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file2
done

