#!/bin/bash
DIR=/home/users/a/avalosma/scratch/12_TRIPLETS/significants
analysis_file=9.1_sizes.txt
echo ' ' > $analysis_file


echo 'signif 1000 perm and FDR' >> $analysis_file
for file in $DIR/*"significant"*.txt ; do
        f=$(echo $file | rev | cut -d "/" -f1 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file
done;

