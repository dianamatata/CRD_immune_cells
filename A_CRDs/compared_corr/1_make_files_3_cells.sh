# for each chromosome, we take the pair of peaks (same for all cells), the corr coeff and the p value for both cell types.
# we add the sample size

out_folder="/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/files"

#sizes

Nneut=165
Nmono=160
Ntcell=94

# for ALL 3 cell types in one file

for chr in {1..22}; do
  
  print $chr
  
  file_mono="/home/users/a/avalosma/scratch/8_PEAKS/hist_mono.corr.chr${chr}.txt"
  file_neut="/home/users/a/avalosma/scratch/8_PEAKS/hist_neut.corr.chr${chr}.txt"
  file_tcell="/home/users/a/avalosma/scratch/8_PEAKS/hist_tcell.corr.chr${chr}.txt"
 
  mono_neut_tcell="$out_folder/mono_neut_tcell.chr${chr}" 
  paste -d" " $file_mono $file_neut $file_tcell | cut -d " " -f1,2,3,4,7,8,15,16,23,24 > ${mono_neut_tcell}.txt


done

# add cis filter!
