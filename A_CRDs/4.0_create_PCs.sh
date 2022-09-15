#!/bin/bash

#DATADIR0=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS
DATADIR1=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells
OUTFOLDER=/home/users/a/avalosma/scratch/4_CRD_residualized
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT $OUTFOLDER/PCs

for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
	BED=$DATADIR1/${cell_type}/qtltools_quantification_filtered.gene.rpkm.bed.gz
	echo "performing QTLtools PCA"
	PCA_FILE=$OUTFOLDER/quantifPCA_${cell_type}
        QTLtools pca --bed $BED --scale --center --out $PCA_FILE

	echo "create PC covariate files"
	for PC in {0..50..2}; do
		cmd="head -$((PC +1)) $PCA_FILE.pca > $OUTFOLDER/PCs/PC_${cell_type}_$PC.txt"
		echo $cmd
		eval $cmd
	done

done

# in DATADIR0 they have the cov sex added at the end

for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
        for PC in {0..50..2}; do
		cat $OUTFOLDER/PCs/${cell_type}_sex_cov.txt >> $OUTFOLDER/PCs/PC_${cell_type}_$PC.txt
	done
done

