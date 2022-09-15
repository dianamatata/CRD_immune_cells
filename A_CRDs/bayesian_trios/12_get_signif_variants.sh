#!/bin/bash
DIR1=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/bayesian_trios
DIR2=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_signif
signif_variants=$DIR2/signif_variants_all.txt

echo ' ' > $signif_variants

for file in $DIR2/*triplet.txt; do
	echo $file
	cat $file| cut -d';' -f1 >> $signif_variants
done

for file in $DIR2/*triplet.txt; do
	f=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | cut -d'_' -f1-2)
	cat $file| cut -d';' -f1 > $DIR2/signif_variant_${f}.txt
done


DIR3=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_all

not_signif_variants=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/bayesian_trios/not_signif_variants

echo ' ' > $not_signif_variants

for file in $DIR3/*; do
        echo $file
        cat $file| cut -d';' -f1 >> $not_signif_variants
done
