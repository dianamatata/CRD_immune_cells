#!/bin/bash

DIR=/home/users/a/avalosma/scratch/5_CRDgene/merged_TH

analysis_file=6.7_signif_analysis_TH.txt
echo ' ' > $analysis_file

for file in $DIR/*.txt.gz ; do
        f=$(echo $file | rev | cut -d "/" -f1 | rev)
        length=$(zcat $file | wc -l)
        echo $length $f >> $analysis_file
done


