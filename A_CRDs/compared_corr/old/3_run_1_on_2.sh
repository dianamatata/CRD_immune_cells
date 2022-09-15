
script_folder=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr
data_folder=/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/files

line=$(cat $data_folder/mono_neut_tcell.chr1.txt | head -1)

Nneut=165
Nmono=160
Ntcell=94

for filename in $data_folder/mono_neut_tcell*.txt; do
  file=$(echo $filename | rev | cut -d'/' -f1 | rev)
  echo $file
  chr=$(echo $file | cut -d'r' -f2 | cut -d'.' -f1)
  
  touch $data_folder/finals/mono_neut_${chr}.txt
  touch $data_folder/finals/mono_tcell_${chr}.txt
  touch $data_folder/finals/neut_tcell_${chr}.txt
  
  while read -r line; do
    peaki=$(echo $line | cut -d' ' -f1)
    peakj=$(echo $line | cut -d' ' -f2)
    loci=$(echo $line | cut -d' ' -f3)
    locj=$(echo $line | cut -d' ' -f4)
    mono_rho=$(echo $line | cut -d' ' -f5)
    mono_pval=$(echo $line | cut -d' ' -f6)
    neut_rho=$(echo $line | cut -d' ' -f7)
    neut_pval=$(echo $line | cut -d' ' -f8)
    tcell_rho=$(echo $line | cut -d' ' -f9)
    tcell_pval=$(echo $line | cut -d' ' -f10)
    

    Rscript $script_folder/old/2_compare_corr.R $mono_rho $mono_pval $neut_rho $neut_pval $Nmono $Nneut $peaki $peakj $loci $locj $data_folder/finals/mono_neut_${chr}.txt
    Rscript $script_folder/old/2_compare_corr.R $mono_rho $mono_pval $tcell_rho $tcell_pval $Nmono $Ntcell $peaki $peakj $loci $locj $data_folder/finals/mono_tcell_${chr}.txt
    Rscript $script_folder/old/2_compare_corr.R $neut_rho $neut_pval $tcell_rho $tcell_pval $Nneut $Ntcell $peaki $peakj $loci $locj $data_folder/finals/neut_tcell_${chr}.txt
    
  done <$filename
done



