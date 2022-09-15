
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

	plink_output=$(echo $OUTDIR/plink/hist_$(echo $filename | cut -d'_' -f2-4)_mean.txt.ld)

        for CRD in $(cat $OUTDIR/temp_CRDs.txt); do
		echo $CRD
                zcat $file | grep $CRD | cut -d' ' -f8 | sort | uniq > $OUTDIR/temp_extract_variants.txt
                count_var_LD=0

		for variant in $(cat $OUTDIR/temp_extract_variants.txt); do
			echo $variant
			cat $plink_output | grep $(echo $variant) | awk '{print $3" "$7}' \
			| tr ' ' '\n' | grep -v $variant > $OUTDIR/temp_matching_variants.txt
				
                        if [ "$(cat $OUTDIR/temp_matching_variants.txt |wc -l)" -ne 0 ]
			then
				vars_in_LD=$(comm -12 <(sort $OUTDIR/temp_extract_variants.txt) <(sort $OUTDIR/temp_matching_variants.txt)) 
				n_vars_in_LD=$(echo $vars_in_LD | grep 'rs' | wc -l)
		
				if [ "$n_vars_in_LD" -ne 0 ]
				then
	                                count_var_LD=$(($count_var_LD+$n_vars_in_LD))
					grep -v $vars_in_LD "$OUTDIR/temp_extract_variants.txt" | grep -v $variant > $OUTDIR/temp
					cat $OUTDIR/temp > $OUTDIR/temp_extract_variants.txt
				fi
			fi
		done

                echo "$CRD $(cat $OUTDIR/temp_extract_variants.txt | wc -l) $count_var_LD" \
		>> $summary_file
	done
done

<< 'MULTILINE-COMMENT'

                        let "n_vars_in_LD=n_vars_in_LD+1" # add +1 for variant
CRD=19_internal_7871
# for each of the variants linked to the CRD,
# extract the variants in LD in the plink file
#remove the vars in LD from pool
# save elements to summary file

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
