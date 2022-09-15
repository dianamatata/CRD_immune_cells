# original version when we need to re,ove the first line of the files, and not here but for file_TF too 

peak_dir="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/peaks_threshold"

data_annot="/home/users/a/avalosma/scratch/Annotations"

file_TF=$data_annot/TFBS_sorted.txt

chr=$1
threshold=$2
file_peak=$peak_dir/peaks_chr_${chr}_threshold_${threshold}.txt
echo peaks_chr_${chr}_threshold_${threshold}.txt

VAR1=$(cat  $file_peak | head -1 | cut -f1)
if [ "$VAR1" = "chr" ]
then
  echo "equal"
  tail -n +2 "$file_peak" > "$file_peak.tmp" && mv "$file_peak.tmp" "$file_peak"
fi

sortBed -i $file_peak > ${file_peak:0:-4}_sorted.txt
rm $file_peak

bedtools intersect -wa -wb \
-a ${file_peak:0:-4}_sorted.txt \
-b $file_TF \
-sorted > $peak_dir/intersect_chr_${chr}_threshold_${threshold}.txt

echo "finished"
