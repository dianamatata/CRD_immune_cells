#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR_0=/home/users/a/avalosma/scratch/0_CRD
DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT $OUT_FOLDER/quantify

for data_type in  'methyl' 'hist' ; do
	for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
        	for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                	data_quantif=$DATADIR_0/quantif_M_${data_type}_${cell_type_quantifM}.bed.gz
                	data_CRD=$DATADIR/${data_type}_${cell_type_CRD}_ALL.CRDs.txt.gz
                	echo "$cell_type_quantifM $cell_type_CRD"
			for c in $(seq 1 22); do
				data_CRD_chr=$DATADIR/${data_type}_${cell_type_CRD}.chr$c\.module.txt.gz
				name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
				LO=$OUT_FOLDER/quantify/$name\.chr$c
	 
				cmd1="$CLOMICs quantify --bed $data_quantif --region $c \
				--tree $data_CRD_chr $data_CRD --out $LO.pc1.txt.gz --pca 1 --normal"
				cmd2="$CLOMICs quantify --bed $data_quantif --region $c \
				--tree $data_CRD_chr $data_CRD --out $LO.mean.txt.gz --mean --normal"
				cmd3="$CLOMICs quantify --bed $data_quantif --region $c \
				--tree $data_CRD_chr $data_CRD --out $LO.loom.txt.gz --loo 0 --normal"
	 
				eval $cmd1; eval $cmd2; eval $cmd3
				JOB=${cell_type}_chr${c}_quantify
				# wsbatch -J $JOB\.job --partition=mono-EL7 --time=00:10:00 \
				# -o $OUT_FOLDER/OUT/$JOB.out -e $OUT_FOLDER/OUT/$JOB.err --wrap="$cmd && $cmd2 && $cmd3"
                	done
	        done
	done
done

