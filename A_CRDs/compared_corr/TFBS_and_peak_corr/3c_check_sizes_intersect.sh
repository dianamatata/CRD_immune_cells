peak_dir="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/peaks_threshold"

for chr in {1..22}
do
  for threshold in 10 50 100 150
  do
    size=$(cat $peak_dir/intersect_chr_${chr}_threshold_${threshold}.txt | wc -l)
    echo "intersect_chr_${chr}_threshold_${threshold}.txt   $size "
  done
done

