#!/bin/bash
DIR_script=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs
INDIR=/home/users/a/avalosma/scratch/7_CRD_Trans/merged
OUTDIR=/home/users/a/avalosma/scratch/7_CRD_Trans/significants
mkdir -p $OUTDIR

for file in $INDIR/*.gz; do
	file=$(echo $file | rev | cut -d "/" -f1 | rev)
	file_no_extension=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '.' -f1)
	Rscript $DIR_script/trans/qtltools_runFDR_ftrans.R $INDIR/$file 0.05 $OUTDIR/FDR_0.05_$file_no_extension
        Rscript $DIR_script/trans/qtltools_runFDR_ftrans.R $INDIR/$file 0.01 $OUTDIR/FDR_0.01_$file_no_extension
done

# need to check which column is taken into account in this file
# also the qvalue can be taken directly in a R script
