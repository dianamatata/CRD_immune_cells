#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR=/home/users/a/avalosma/scratch/2_CRD
DATA_TREE=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/8_CRD_PEAKS
mkdir -p $OUT_FOLDER $OUT_FOLDER/chr $OUT_FOLDER/ALL

# LM is only the list of CRD names, so look at 7.1 first
# module does not apply here. aCRD or sCRD: it is still the same CRD ID that we extract

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' ; do
			LM=$DATADIR/CRD_names/${data_type}_${cell_type}.ALLchr.${module}.txt
			name=${data_type}_${cell_type}
			for c in $(seq 1 22); do
				LT=$DATA_TREE/${data_type}_${cell_type}.chr$c\.tree.txt.gz
			        LO=$OUT_FOLDER/chr/${name}.chr$c\.peak
        		        cmd="$CLOMICs gpeak --tree $LT --mod $LM --out ${LO}.txt.gz"
                		eval $cmd
				# JOB=chr${c}_peak_$name
				# wsbatch -J $JOB\.job --partition=mono-shared-EL7 --time=00:00:15 --wrap="$cmd"
			done
		done
	done
done


# merge
for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean'  ; do
		name=${data_type}_${cell_type}
			for c in $(seq 1 22); do
	                	LO=$OUT_FOLDER/chr/${name}.chr$c\.peak
        	        	LO_ALL=$OUT_FOLDER/ALL/${name}.ALLchr.peaks.txt
                		zcat $LO\.txt.gz >> $LO_ALL
			done
                done
        done
done


