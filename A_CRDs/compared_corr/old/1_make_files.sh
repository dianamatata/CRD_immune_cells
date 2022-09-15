# for each chromosome, we take the pair of peaks ( same for all cells ), the corr coeff and the p value for both cell types.
# we add the sample size

out_folder="/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/files"

Nneut=165
Nmono=160
Ntcell=94



for chr in {1..22}; do
  
  print $chr
  
  file_mono="/home/users/a/avalosma/scratch/8_PEAKS/hist_mono.corr.chr${chr}.txt"
  file_neut="/home/users/a/avalosma/scratch/8_PEAKS/hist_neut.corr.chr${chr}.txt"
  file_tcell="/home/users/a/avalosma/scratch/8_PEAKS/hist_tcell.corr.chr${chr}.txt"
 
  mono_neut_tcell="$out_folder/mono_neut_tcell.chr${chr}" 
  paste -d" " $file_mono $file_neut $file_tcell | cut -d " " -f1,2,3,4,7,8,15,16,23,24 > ${mono_neut_tcell}.txt

  # create input files
  mono_neut="$out_folder/mono_neut.chr${chr}"
  paste -d" " $file_mono $file_neut | cut -d " " -f7,8,15,16 > ${mono_neut}.txt
  
  mono_tcell="$out_folder/mono_tcell.chr${chr}"
  paste -d" " $file_mono $file_tcell | cut -d " " -f7,8,15,16 > ${mono_tcell}.txt
  
  neut_tcell="$out_folder/neut_tcell.chr${chr}"
  paste -d" " $file_neut $file_tcell | cut -d " " -f7,8,15,16 > ${neut_tcell}.txt
  
  echo "finished"

done


#  paste <(cat $file_mono  | awk '{print $7 "\t" $8}') <(cat $file_neut | awk '{print $7 "\t" $8}') | sed 's/$/\t160\t165/' > ${mono_neut}.txt
#  paste <(cat $file_mono  | awk '{print $7 "\t" $8}') <(cat $file_tcell | awk '{print $7 "\t" $8}') | sed 's/$/\t160\t94/' > ${mono_tcell}.txt
#  paste <(cat $file_neut  | awk '{print $7 "\t" $8}') <(cat $file_tcell | awk '{print $7 "\t" $8}') | sed 's/$/\t165\t94/' > ${neut_tcell}.txt


