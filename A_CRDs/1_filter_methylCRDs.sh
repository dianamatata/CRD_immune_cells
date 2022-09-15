#!/bin/bash
# filter the CRDs

DATADIR=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs
OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD
summary_file=$OUT_FOLDER/summary_methylCRDs
echo " cell_type post_CRDs threshold pre_CRDs" > $summary_file

for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do

	cmd="Rscript $DATADIR/1b_determine_methylCRD_threshold.R '$cell_type' $OUT_FOLDER"
	eval $cmd
        threshold=$(cat $OUT_FOLDER/threshold_methylCRD_$cell_type)
	threshold_int=$(echo $threshold | awk '{printf("%d\n",$1 + 0.5)}')
        echo -e $cell_type ' \t ' $threshold ' \t ' $threshold_int

	for c in $(seq 1 22); do
                LO=$OUT_FOLDER/${cell_type}.chr$c
                cmd="zcat ${LO}.module.txt.gz | awk '{ if (\$30 == 1 && \$19 > $threshold_int) print \$4}'"
                eval $cmd
	done | gzip -c > $OUT_FOLDER/$cell_type\_ALL.CRDs.txt.gz
	
	echo " ${cell_type} $(zcat $OUT_FOLDER/$cell_type\_ALL.CRDs.txt.gz | wc -l) $threshold_int $(cat $OUT_FOLDER/${cell_type}_dist.ALLchr.bed | wc -l)" >> $summary_file
done
