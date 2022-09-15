#!/bin/bash
permuts=1000
DIR_script=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/bayesian_trios
INDIR=/home/users/a/avalosma/scratch/12_TRIPLETS/merged_$permuts
OUTDIR=/home/users/a/avalosma/scratch/12_TRIPLETS/significants
mkdir -p $OUTDIR

echo ' ' > $OUTDIR/files.txt

for file in $INDIR/*.txt.gz; do
	file=$(echo $file | rev | cut -d "/" -f1 | rev)
	filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '.' -f1)
	echo $file >> $OUTDIR/files.txt
	Rscript $DIR_script/2.1_qtltools_runFDR_cis.R $INDIR/$file 0.05 $OUTDIR/FDR_0.05_$filename
        Rscript $DIR_script/2.1_qtltools_runFDR_cis.R $INDIR/$file 0.01 $OUTDIR/FDR_0.01_$filename
done
