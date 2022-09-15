#! /bin/bash

outpath='/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS'
cd $outpath

for compare_file in 'COR_mono_vs_neut_chr' 'COR_mono_vs_tcell_chr' 'COR_neut_vs_tcell_chr'; do
	echo $compare_file
	for chr in {1..22}; do
		ls $outpath/${compare_file}${chr}_section*.txt | xargs -n 1 tail -n +2 | > $outpath/merged/${compare_file}${chr}_ALL.txt
		cat $outpath/merged/${compare_file}${chr}_ALL.txt  | cut -f2- | sort -k2 -n > $outpath/merged/${compare_file}${chr}_sorted.txt
	done
done

