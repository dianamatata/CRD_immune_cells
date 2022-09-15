#!/bin/bash

DATADIR0=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS
DATADIR1=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells
OUTFOLDER=/home/users/a/avalosma/scratch/4_CRD_residualized
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

# for RNA
for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
	for PC in {0..50..2}; do
		INP=$DATADIR1/${cell_type}/qtltools_quantification_filtered.gene.rpkm.bed.gz
		COV=$DATADIR0/${cell_type}_v3.0/covariates_$PC.txt
		OUT=$OUTFOLDER/${cell_type}_RNA.PC$PC\.bed
		
		cmd="QTLtools correct --bed $INP --cov $COV --normal --out $OUT && bgzip $OUT && tabix -p bed $OUT\.gz"
		eval $cmd
		#JOB=residualize_$PC
		#wsbatch -J $JOB\.job --partition=mono-EL7 --time=00:01:00 -o $OUT_FOLDER/OUT/$JOB.out -e $OUT_FOLDER/OUT/$JOB.err --wrap="$cmd"
	done
done

