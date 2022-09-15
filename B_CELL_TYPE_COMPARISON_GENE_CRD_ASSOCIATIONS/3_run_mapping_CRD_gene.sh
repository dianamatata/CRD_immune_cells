#!/bin/bash
# in the crossed version we want to mix the cell types 
# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
# for the moment no reverse and the TH does not work

PC=10
K=100
DATADIR0=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS
DATADIR=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
DIR_RNA=/home/users/a/avalosma/scratch/4_CRD_residualized
OUTDIR=$DATADIR/mapping_gene_CRDs

mkdir -p $OUTDIR $OUTDIR/mean $OUTDIR/loom
for module in 'mean' 'loom' ; do
        mkdir -p $OUTDIR/${module}/mapping_CRD_gene
done

declare -A rna_file
rna_file[neut]="EGAD00001002675"
rna_file[mono]="EGAD00001002674"
rna_file[tcell]="EGAD00001002671"


for data_type in  'methyl' 'hist' ; do
        for cell_type_quantifM in 'neut'  'mono' 'tcell' ; do
                for cell_type_CRD in 'neut'  'mono' 'tcell' ; do
                        name=${data_type}_${cell_type_quantifM}_vs_${cell_type_CRD}
                        for module in 'mean' 'loom' ; do
        			MOD=$DATADIR/quantify_ALL/${name}.ALLchr.${module}.txt.gz
        	                RNA=$DIR_RNA/${rna_file[${cell_type_quantifM}]}_RNA.PC$PC\.bed.gz
				for k in $(seq 1 $K); do
				        OUT1=$OUTDIR/${module}/mapping_CRD_gene/${name}_CRD_gene_chunk$k
	                                cmd1="QTLtools cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out ${OUT1}.txt"
		                        eval $cmd1
				done
			done
		done
	done
done

