#!/bin/bash

DATADIR=/home/users/a/avalosma/scratch/4_CRD_residualized
DATACRD=/home/users/a/avalosma/scratch/2_CRD
DIR_2col=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1000/2_columns
CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
CLOMICs=/srv/beegfs/scratch/groups/funpopgen/Tools/clomics1/bin/clomics
DIR_2col=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1000/2_columns
DIR_2col=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1/2_columns
OUT=/home/users/a/avalosma/scratch/12_TRIPLETS/PC1
OUT=/home/users/a/avalosma/scratch/12_TRIPLETS/PC1_all
mkdir -p $OUT $OUT/quantify


declare -A rna_file
rna_file[neut]="EGAD00001002675"
rna_file[mono]="EGAD00001002674"
rna_file[tcell]="EGAD00001002671"
PC=10

module='mean'
for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		RNA=$DATADIR/${rna_file[${cell_type}]}_RNA.PC$PC\.bed.gz
		CRD=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
		OUTPUT=$OUT/quantify/${data_type}_${cell_type}
		COL2=$DIR_2col/${data_type}_${cell_type}_mean_mapping_CRD_gene_ALL.txt
		cmd="$CLOMICs quantify --bed $RNA $CRD --grp $COL2 --pca 1 --normal --out $OUTPUT.pc1.txt"
		eval $cmd
        done
done



# << 'MULTILINE-COMMENT'
for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                OUTPUT=$OUT/quantify/${data_type}_${cell_type}
		awk 'NR == 1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' $OUTPUT.pc1.txt \
		 | bgzip -c > $OUTPUT.pc1.sorted.txt.gz
		tabix -p bed $OUTPUT.pc1.sorted.txt.gz
	done
done

# MULTILINE-COMMENT

