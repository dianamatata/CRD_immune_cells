peak_dir="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/peaks_threshold"
data_annot="/home/users/a/avalosma/scratch/Annotations"

file_TF=$data_annot/TFBS_sorted.txt

chr=$1
threshold=$2

file_peak=$peak_dir/peaks_chr_${chr}_threshold_${threshold}.txt
echo peaks_chr_${chr}_threshold_${threshold}.txt

# sortBed -i $file_peak > ${file_peak:0:-4}_sorted.txt

bedtools intersect -wa -wb \
-a ${file_peak:0:-4}_sorted.txt \
-b $file_TF \
-sorted > $peak_dir/intersect_chr_${chr}_threshold_${threshold}.txt

echo "finished"
