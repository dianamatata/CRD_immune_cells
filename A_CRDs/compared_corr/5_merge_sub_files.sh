folder='/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/files'

for chr in {1..22}
do
        echo $chr
        ls $folder/sub_files2/CORR_mono_neut_tcell_chr_${chr}_section*.txt | xargs -n 1 tail -n +2 > $folder/merged/CORR_mono_neut_tcell_chr_${chr}_ALL.txt
        cat $folder/merged/CORR_mono_neut_tcell_chr_${chr}_ALL.txt | sort -n -k1,1 -k2,2 > $folder/merged/CORR_mono_neut_tcell_chr_${chr}_ALL_sorted.txt
#        rm $folder/merged/CORR_mono_neut_tcell_chr_${chr}_ALL.txt
done

rm *ALL.txt
