FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
OUTFOLDER="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TF_arrays"
for chr in {1..22}
do
	for pair in "mono_neut" "mono_tcell" "neut_tcell" "neut_mono" "tcell_mono" "tcell_neut"
	do
		file=$OUTFOLDER/${pair}_${chr}_array_specific.txt
		if [ -f "$file" ]; then
		    echo "${pair}_${chr} exists."
		else 
		    echo "computing ${pair}_${chr} "
			cmd="Rscript $FOLDER/TFBS_and_peak_corr/7_compute_arrays_for_heatmap.R $chr $pair"
			echo $cmd
			wsbatch -J heatmap_${chr}_${pair}.job \
			--partition=public-cpu \
			--time=24:00:00 \
			--mem=10000 \ 
			-o $FOLDER/log/heatmap_${chr}_${pair}.out \ 
			-e $FOLDER/log/heatmap_${chr}_${pair}.err \
			--wrap="$cmd"
		fi
	done
done

