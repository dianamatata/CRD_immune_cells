#!/bin/bash
DIR_script=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs
INDIR=/home/users/a/avalosma/scratch/5_CRDgene/merged
OUTDIR=/home/users/a/avalosma/scratch/5_CRDgene/significants
mkdir -p $OUTDIR

# compute qval from col 20: ~/IMMUNE_CELLS/diana_scripts/A_CRDs
# check man qtltools cis
for file in $INDIR/*mapping_CRD_gene_ALL.txt.gz; do
	file=$(echo $file | rev | cut -d "/" -f1 | rev)
	filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '.' -f1)
	Rscript $DIR_script/6.1_qtltools_runFDR_cis.R $INDIR/$file 0.05 $OUTDIR/FDR_0.05_$filename
        Rscript $DIR_script/6.1_qtltools_runFDR_cis.R $INDIR/$file 0.01 $OUTDIR/FDR_0.01_$filename
done
