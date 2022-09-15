#!/bin/bash
DIR_script=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs
INDIR=/home/users/a/avalosma/scratch/5_CRDgene/merged_TH
OUTDIR=/home/users/a/avalosma/scratch/5_CRDgene/significants_TH
mkdir -p $OUTDIR

for file in $INDIR/*.gz; do
	file=$(echo $file | rev | cut -d "/" -f1 | rev)
	filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '.' -f1)
	Rscript $DIR_script/6.1_qtltools_runFDR_cis.R $INDIR/$file 0.05 $OUTDIR/FDR_0.05_$filename
        Rscript $DIR_script/6.1_qtltools_runFDR_cis.R $INDIR/$file 0.01 $OUTDIR/FDR_0.01_$filename
done
