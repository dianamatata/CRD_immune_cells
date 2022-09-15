#!/bin/bash

#run QTL cis permutations for all the number of PCs and then select the number of PCs which maximizes aCRD and sCRD-QTLs discovery

DATACRD=/home/users/a/avalosma/scratch/2_CRD
DATADIR=/home/users/a/avalosma/scratch/4_CRD_residualized
mkdir -p $DATADIR/mapping_aCRDs_many_PCs $DATADIR/mapping_sCRDs_many_PCs

OUT_aCRDs=$DATADIR/mapping_aCRDs_many_PCs
OUT_sCRDs=$DATADIR/mapping_sCRDs_many_PCs

# look for equivalent names of cell types from chip seq and methyl marks to RNA data
declare -A cell_types
cell_types[EGAD00001002670]="EGAD00001002675"
cell_types[EGAD00001002672]="EGAD00001002674"
cell_types[EGAD00001002673]="EGAD00001002671"
cell_types[methyl_neut]="EGAD00001002675"
cell_types[methyl_mono]="EGAD00001002674"
cell_types[methyl_tcell]="EGAD00001002671"

permutations=10000
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
        aCRD=$DATACRD/${cell_type}.ALLchr.mean.txt.gz
	sCRD=$DATACRD/${cell_type}.ALLchr.loom.txt.gz
	for PC in {0..50..2}; do
		RNA=$DATADIR/${cell_types[${cell_type}]}_RNA.PC$PC\.bed.gz        
		COV=$DATADIR/PCs/PC_${cell_types[${cell_type}]}_$PC.txt
               	cmd_aCRD_cis="QTLtools cis --vcf $aCRD --bed $RNA --cov $COV --permute $permutations --out $OUT_aCRDs/permut_${cell_type}_$PC.txt"
                cmd_sCRD_cis="QTLtools cis --vcf $sCRD --bed $RNA --cov $COV --permute $permutations --out $OUT_sCRDs/permut_${cell_type}_$PC.txt"
		echo $cmd_aCRD_cis
                JOB=${cell_type}_${PC}_aCRD_$permutations
                wsbatch -J $JOB\.job --partition=mono-EL7 --time=10:00:00 -o $DATADIR/OUT_many_PCs/$JOB.out -e $DATADIR/OUT_many_PCs/$JOB.err --wrap="$cmd_aCRD_cis"
                JOB=${cell_type}_${PC}_sCRD_$permutations
                wsbatch -J $JOB\.job --partition=mono-EL7 --time=10:00:00 -o $DATADIR/OUT_many_PCs/$JOB.out -e $DATADIR/OUT_many_PCs/$JOB.err --wrap="$cmd_sCRD_cis"
        done
done

# compute eQTLs, methQTLs, histQTLs, aCRD and sCRD-QTLs. 
# need to define the optimal number of PC for that?

#               time 10:00 job canceled due to time limit, but methyl_neut, PC50 takes 500sec ~ 8min
#               eval $cmd_aCRD_cis
#               eval $cmd_sCRD_cis
