#!/bin/bash

DIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/conditional_merged
OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/LD

# for all the CRD-QTLs after the conditional pass
# extract for each CRD the variants associated
# are these variants in LD?

for file in $DIR/*.gz; do
        zcat $file | cut -d' ' -f1 | uniq > $OUTDIR/temp_CRDs.txt
        filename=$( echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
	echo $filename
        summary_file=$OUTDIR/${filename}_summary.txt
        echo "CRD_name	N_variants_associated	N_variants_in_LD" > $summary_file
        plink_output=$(echo $OUTDIR/plink/$(echo $filename | cut -d'_' -f1-4)_mean.txt.ld)
       
	 for CRD in $(cat $OUTDIR/temp_CRDs.txt); do
	        zcat $file | grep $CRD | cut -d' ' -f8 > $OUTDIR/temp_extract_variants.txt
	
		# for each of the variants linked to the CRD,
		# extract the variants in LD in the plink file
                for variant in $(cat $OUTDIR/temp_extract_variants.txt); do
			cat $plink_output | grep $(echo $variant) | awk '{print $3" "$7}' \
			| tr ' ' '\n' | grep -v $variant > $OUTDIR/temp_matching_variants.txt
		
			n_vars_in_LD=$(comm -12 <(sort $OUTDIR/temp_extract_variants.txt) <(sort $OUTDIR/temp_matching_variants.txt) | wc -l)

			# save elements to summary file
                	echo "$CRD $(cat $OUTDIR/temp_extract_variants.txt | wc -l) $n_vars_in_LD" >> $summary_file
        	done
	done
done


# for each of the variants linked to the CRD,
# extract the variants in LD in the plink file

find_variants_in_LD_with_CRD() {
	echo $CRD
	for variant in $(cat $OUTDIR/temp_extract_variants.txt); do
		cat $plink_output | grep $(echo $variant) | awk '{print $3" "$7}' \
		| tr ' ' '\n' | grep -v $variant > $OUTDIR/temp_matching_variants.txt

		n_vars_in_LD=$(comm -12 <(sort $OUTDIR/temp_extract_variants.txt) <(sort $OUTDIR/temp_matching_variants.txt) | wc -l)

		# save elements to summary file
		echo "$CRD $(cat $OUTDIR/temp_extract_variants.txt | wc -l) $n_vars_in_LD" >> $summary_file
	done

}

<< 'MULTILINE-COMMENT'
                       if [ "$n_vars_in_LD" -ne 0 ]; then echo $(comm -12 <(sort $OUTDIR/temp_extract_variants.txt) <(sort $OUTDIR/temp_matching_variants.txt)); fi
#debug
file=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/conditional_merged/hist_mono_vs_neut_mean.txt.gz
CRD=12_internal_14191

                # extract variants linked to CRD that are in the plink file
                echo ' ' > $OUTDIR/grep_variants_CRD_LD.txt
                for variant in $(cat $OUTDIR/temp_extract_variants.txt); do
                        grep "$variant" $plink_output | awk '{print $3" "$7" "$9}' >> $OUTDIR/grep_variants_CRD_LD.txt
                done
MULTILINE-COMMENT
