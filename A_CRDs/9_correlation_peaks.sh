#!/bin/bash
# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/bin/run_clomics_corr_v3.0.py

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR_0=/home/users/a/avalosma/scratch/0_CRD
DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/8_PEAKS

mkdir -p $OUT_FOLDER

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                LI=$DATADIR_0/quantif_M_${data_type}_${cell_type}.bed.gz
                ### Residualise
		# QTLtools correct --bed $bedfile --cov $covariate --norm --out $bed_residuals

		### Run clomics correlations
		for c in $(seq 1 22); do
                        LO=$OUT_FOLDER/${data_type}_${cell_type}.corr.chr$c
			cmd="$CLOMICs corr --bed $LI --out $LO.txt --region $c"
			wsbatch -J peak_${data_type}_${cell_type}_${c}.job --partition=mono-shared-EL7 --time=04:00:00 --wrap="$cmd"
		done
	done
done


