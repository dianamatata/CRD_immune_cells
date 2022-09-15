#!/bin/bash

#QTLTOOLS=/home/users/a/avalosma/bin/qtltools
IN=/home/users/a/avalosma/scratch/Blueprint/ref_files/gencode.v15.annotation.gtf.gz
IN_BAM=/home/users/a/avalosma/scratch/Blueprint/RNA_seq
OUTFOLDER=/home/users/a/avalosma/scratch/RNA_seq
mkdir -p $OUTPUT_FOLDER $OUTPUT_FOLDER/OUT

#for cell_type in 'EGAD00001002675'  'EGAD00001002674' 'EGAD00001002671' ; do 
cell_type='EGAD00001002675'

cmd="QTLtools quan --gtf $IN --bam $IN_BAM/$cell_type/*.bam --out qtltools_quantification_$cell_type --filter-mapping-quality 255 --rpkm "

echo $cmd

JOB=qtltools_quant_rpkm

sbatch -J $JOB\.job --partition=mono-EL7 --time=20:00:00 -o $OUT_FOLDER/$JOB.out -e $OUT_FOLDER/$JOB.err --wrap="$cmd"
