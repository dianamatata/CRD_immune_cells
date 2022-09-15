#!/bin/bash
# filter the CRDs
# cond2 mean the thresholds set by Guillaume, which are a bit different...

DATADIR=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs
OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD
mkdir -p $OUTFOLDER/cond2
summary_file=$OUT_FOLDER/cond2/summary_methylCRDs
echo " cell_type post_CRDs threshold pre_CRDs" > $summary_file

declare -A thresholds
thresholds[methyl_neut]=3800
thresholds[methyl_mono]=3300
thresholds[methyl_tcell]=4300

for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
        threshold_int=$thresholds[${cell_type}]
        echo -e $cell_type ' \t ' $threshold_int
        for c in $(seq 1 22); do
	        LO=$OUT_FOLDER/${cell_type}.chr$c
        	cmd="zcat ${LO}.module.txt.gz | awk '{ if (\$30 == 1 && \$19 > $threshold_int) print \$4}'"
		eval $cmd
        done | gzip -c > $OUT_FOLDER/cond2/$cell_type\_ALL.CRDs.txt.gz
	echo " ${cell_type} $(zcat $OUT_FOLDER/cond2/$cell_type\_ALL.CRDs.txt.gz | wc -l) $threshold_int $(cat $OUT_FOLDER/${cell_type}_dist.ALLchr.bed | wc -l)" >> $summary_file
done
