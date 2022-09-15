FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"

for chr in {1..22}
do
  for sub in {1..10}
  do
	cmd="Rscript $FOLDER/3_compare_corr_coeff_btw_peaks_v2.R $chr $sub"
	wsbatch -J Peak_corr_sub_${chr}.job \
	--partition=shared-bigmem  \
	--time=12:00:00 \
	--mem=100000  \
	-o $FOLDER/log/Peak_corr_sub_${chr}.out  \
	-e $FOLDER/log/Peak_corr_sub_${chr}.err  \
	--wrap="$cmd"
  done
done
