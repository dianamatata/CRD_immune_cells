#!/bin/bash

### merge all names of shared CRDs in one list

DIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/CRD_sharing

outfile=/home/users/a/avalosma/scratch/10_CRD_QTLs/CRD_sharing/shared_CRDs_temp.txt
outfile2=/home/users/a/avalosma/scratch/10_CRD_QTLs/CRD_sharing/shared_CRDs_ALL.txt

for file in $DIR/*sharedCRDs.txt; do
	cat $file | sed 's/ /\n/g' >> $outfile
done

cat $outfile | sort | uniq > $outfile2
rm $outfile


### filter all triplets in /home/users/a/avalosma/scratch/12_TRIPLETS/not_signif/merged_nominal1/triplets/ that have a shared CRD name and a significant variant

DIR_TRIPLETS=/home/users/a/avalosma/scratch/12_TRIPLETS/not_signif/merged_nominal1/triplets
signif_variants=/home/users/a/avalosma/scratch/12_TRIPLETS/vcf/all_signif_variants.txt

for file in $DIR_TRIPLETS/*ALL.txt; do
	echo $file
	f=$(echo $file| rev | cut -d'/' -f1 | rev)
	for line in $(cat $file); do
		echo $line
	done
done

#	memory exhausted
#	grep -w -F -f $outfile2 $file > $DIR_TRIPLETS/filtered_triplets_common_CRDs/$f
#	grep -w -F -f $signif_variants $(grep -w -F -f $outfile2 $file) > $DIR_TRIPLETS/filtered_triplets_common_CRDs/$f
done

: <<'END'

for file in $DIR_TRIPLETS/filtered_triplets_common_CRDs/*ALL.txt; do
	f=$(echo $file| rev | cut -d'/' -f1 | rev)
	grep -w -F -f $signif_variants $file > $DIR_TRIPLETS/filtered_triplets_common_CRDs_variants/$f
done
END
