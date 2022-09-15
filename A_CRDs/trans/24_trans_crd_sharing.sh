#!/bin/bash
# invalid compressed data--format violated
# the trans files are not right
# try with trans crd files from guillaume first
# them understand what is wrong with these

### FOLDERS
OUTDIR=/home/users/a/avalosma/scratch/7_CRD_Trans
baobab_trans_dir=/home/users/a/avalosma/scratch/7_CRD_Trans/merged
mac_trans_dir=/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant

### FILES
shared_crds=/Users/dianaavalos/Programming/CRDs_nov20/debug/analysis_files/mono_vs_neut_sharedCRDs.txt
trans_crd_mono=${mac_trans_dir}/hist_mono_mean_trans.significant_0.01.txt
trans_crd_neut=${mac_trans_dir}/hist_neut_mean_trans.significant_0.01.txt

### GETS PAIRS ONLY
cat $trans_crd_mono | cut -d ' ' -f4,8 > $mac_trans_dir/trans_mono_crds_corr.txt
cat $trans_crd_neut | cut -d ' ' -f4,8 > $mac_trans_dir/trans_neut_crds_corr.txt


# query ref, query is mono, ref is neut
# mono vs neut

### EXTRACT SUBSET OF CRDs OVERLAPPING IN EACH CELL TYPE
crds_mono_overlapping=$(cat $shared_crds | cut -d ' ' -f1)
crds_neut_overlapping=$(cat $shared_crds | cut -d ' ' -f2)

### WRITE SUBSET OF CRDs OVERLAPPING IN EACH CELL TYPE
cat $shared_crds | cut -d ' ' -f1 > ${mac_trans_dir}/crds_mono_o.txt
cat $shared_crds | cut -d ' ' -f2 > ${mac_trans_dir}/crds_neut_o.txt

### ONLY PAIRS OF CRDs FROM TRANS CRD FILE
gunzip -c $trans_crd_mono | head | cut -d ' ' -f4,8
gunzip -c $trans_crd_mono | cut -d ' ' -f4,8 > trans_crds_mono.txt 
# uncompress failed , try in baobab

# weird why is the output different?
cat $shared_crds


${mac_trans_dir}/hist_mono_mean_trans.subset.txt
fgrep -w -o -f $(cat $shared_crds | cut -d ' ' -f1) ${mac_trans_dir}/hist_mono_mean_trans.subset.txt > out.txt

# extract words from file1 in file2, different possibilities
$(cat $shared_crds | cut -d ' ' -f1) 
fgrep -w -o -f $crds_mono_overlapping $(gunzip -c $trans_crd_mono) > out.txt
grep -w -Ff  $crds_mono_overlapping $(gunzip -c $trans_crd_mono) > out2.txt
grep -o -f file1 file2
comm -12 list.txt output.txt
comm -12 <(sort list.txt) <(sort output.txt) 
# faster if sorted
grep -wFf list.txt output.txt 
 

$(cat $shared_crds | cut -d ' ' -f1) $(gunzip -c $trans_crd_mono) > out.txt
'20_internal_5117'  $(gunzip -c $trans_crd_mono) > out.txt 



# in the transCRD of neut, we subselect the crds that are common with mono (take neut column), with grep
# replace the neut crds with their correspondant? 
# problem if there are many matches? the many matches are already duplicated in the correspondance shared_crds
# we do the same in mono
# we look how many pairs are in both: grep from a file to another

gunzip -c $trans_crd_mono# we have subsets of shared_crds with the common crds, and 2 crds per line that are associated
# we compare the lines that match


: <<'END'
1 10 121204 10_internal_10433 1 11 222036 11_internal_8165 -0.142 0.0741626
1 10 121204 10_internal_10433 2 11 283161 11_internal_8292 -0.003 0.966251
END


# grep from 1 shared_crds to another
grep -w -Ff shared_crds1 file2
fgrep -w -o -f "words.txt" "text.txt"

: <<'END'
-w     Ne sélectionner que les lignes contenant une correspondance  formant  un mot complet.
-f fichier Lire le motif dans le fichier indiqué
-F     Interprète le motif  comme  une  liste  de  chaînes figées, séparées par des Sauts de Lignes
END

# we have the CRD matching shared_crds / pairs of cell types / data_type
