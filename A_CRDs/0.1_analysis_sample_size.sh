#!/bin/bash

OUT_FOLDER=/home/users/a/avalosma/scratch/0_CRD
ANALYSIS_FILE=$OUT_FOLDER/analysis_sample_size.txt

'' > $ANALYSIS_FILE

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
        sample_count=$(( $(zcat $OUT_FOLDER/${cell_type}_merged_residuals.bed.gz | head -1 | wc -w) -6 ))
        echo -e $cell_type' \t '$sample_count >> $ANALYSIS_FILE
done
for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
	sample_count=$(( $(zcat $OUT_FOLDER/${cell_type}.bed.gz | head -1 | wc -w) -6 ))
        echo -e $cell_type' \t '$sample_count >> $ANALYSIS_FILE
done
