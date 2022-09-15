#! /bin/bash

FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"

cmd="Rscript $FOLDER/compare_corr_coeff.R"
wsbatch -J compute_corr.job --partition=shared-bigmem --time=12:00:00 --mem=100000 -o $FOLDER/log/compute_corr.out -e $FOLDER/log/compute_corr.err --wrap="$cmd"


for chr in {1..22}
do
	cmd="Rscript $FOLDER/compare_corr_coeff_per_chr.R $chr"
	echo $cmd
	wsbatch -J corr_$chr\.job --partition=shared-bigmem --time=12:00:00 --mem=100000 -o $FOLDER/log/corr_$chr\.out -e $FOLDER/log/corr_$chr\.err --wrap="$cmd"
done

cmd2="Rscript $FOLDER/compare_corr_coeffs_btw_CRDs.R"
wsbatch -J corr_CRD.job --partition=shared-bigmem --time=12:00:00 --mem=100000 -o $FOLDER/log/corr_CRD.out -e $FOLDER/log/corr_CRD.err --wrap="$cmd2"


cmd="Rscript $FOLDER/compare_corr_coeff_per_chr_v2.R"
wsbatch -J v2.job --partition=shared-bigmem --time=12:00:00 --mem=100000 -o $FOLDER/log/v2.out -e $FOLDER/log/v2.err --wrap="$cmd"



cmd="Rscript $FOLDER/compare_coeff_per_chr.sh"
wsbatch -J vsh.job --partition=shared-bigmem --time=12:00:00 --mem=100000 -o $FOLDER/log/vsh.out -e $FOLDER/log/vsh.err --wrap="$cmd"


# last version
for chr in {1..22}
do
	cmd="Rscript $FOLDER/compare_corr_coeff_btw_peaks.R $chr"
	wsbatch -J peak_corr_${chr}.job --partition=public-cpu --time=24:00:00 --mem=10000  -o $FOLDER/log/peak_corr_${chr}.out -e $FOLDER/log/peak_corr_${chr}.err --wrap="$cmd"
done

# failed for 17 12 11 10 7 6 5 3 2 1

for chr in 17 12 11 10 7 6 5 3 2 1
do
        cmd="Rscript $FOLDER/compare_corr_coeff_btw_peaks.R $chr"
        wsbatch -J peak_corr_${chr}.job --partition=public-cpu --time=2-24:00:00 --mem=10000  -o $FOLDER/log/peak_corr_${chr}.out -e $FOLDER/log/peak_corr_${chr}.err --wrap="$cmd"
done
