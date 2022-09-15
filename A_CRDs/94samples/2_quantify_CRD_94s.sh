#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR_0=/home/users/a/avalosma/scratch/0_CRD/94_samples
DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/2_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	LI=$DATADIR_0/${cell_type}_merged_residuals_94s.bed.gz
        LM=$DATADIR/${cell_type}_94s_ALL.modules.MOD1.NRE2.txt.gz
	for c in $(seq 1 22); do
		echo "$cell_type $c"
		LT=$DATADIR/${cell_type}_94s.chr$c\.module.txt.gz
		LO=$OUT_FOLDER/${cell_type}_94s.chr$c\

		cmd="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.pc1.txt.gz --pca 1 --normal"
		cmd2="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.mean.txt.gz --mean --normal"
		cmd3="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.loom.txt.gz --loo 0 --normal"

                eval $cmd; eval $cmd2; eval $cmd3
                JOB=$OUT_FOLDER/OUT/${cell_type}_chr${c}_quantify_94s
		#sbatch -J $JOB\.job --partition=mono-EL7 --time=00:01:00 -o $OUT_FOLDER/OUT/$JOB.out -e $OUT_FOLDER/OUT/$JOB.err --wrap="$cmd && $cmd2 && $cmd3" 
	done
done
