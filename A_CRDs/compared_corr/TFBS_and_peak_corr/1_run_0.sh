#!/bin/bash
FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"

cmd="Rscript $FOLDER/TFBS_and_peak_corr/0_create_TFBS_file.R"
echo $cmd
wsbatch -J TF_file \
--partition=shared-bigmem \
--time=12:00:00 \
--mem=100000 \
-o $FOLDER/log/TF_file.out \
-e $FOLDER/log/TF_file.err \
--wrap="$cmd"



cmd2="Rscript $FOLDER/TFBS_and_peak_corr/2_create_peak_files.R"

echo $cmd3
wsbatch -J peak_file \
--partition=shared-cpu \
--time=12:00:00 \
--mem=10000 \
-o $FOLDER/log/peak_file.out \
-e $FOLDER/log/peak_file.err \
--wrap="$cmd2"

cmd3="source $FOLDER/TFBS_and_peak_corr/3_bed_intersect.sh"
echo $cmd3
wsbatch -J intersect2 \
--partition=shared-cpu \
--time=12:00:00 \
--mem=10000 \
-o $FOLDER/log/intersect2.out \
-e $FOLDER/log/intersect2.err \
--wrap="$cmd3"


for chr in {1..22}
do
  for threshold in 10 50 100 150
  do
    print $chr $threshold
    cmd4="Rscript $FOLDER/TFBS_and_peak_corr/4_create_freq_TF_per_bin.R $chr $threshold"
    wsbatch -J bins_${chr}_${threshold} \
    --partition=shared-cpu \
    --time=12:00:00 \
    --mem=10000 \
    -o $FOLDER/log/bins_${chr}_${threshold}.out \
    -e $FOLDER/log/bins_${chr}_${threshold}.err \
    --wrap="$cmd4"
  done
done


