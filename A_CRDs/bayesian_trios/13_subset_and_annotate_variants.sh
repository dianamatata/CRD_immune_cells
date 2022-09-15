#!/bin/bash

VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
BCFTOOLS=/srv/beegfs/scratch/groups/funpopgen/Tools/bcftools-1.10.2/bcftools
DIR=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_signif
DIR_VCF=/home/users/a/avalosma/scratch/12_TRIPLETS/vcf
DIR_VCF_A=$DIR_VCF/annotated
mkdir -p $DIR_VCF $DIR_VCF_A

VCFfile_signif_annotated=$DIR_VCF/VCFfile_signif_annotated.vcf.gz
VCFfile_signif_annotated2=$DIR_VCF/VCFfile_signif_annotated2.vcf.gz

#### create vcf files for each file of signif variant per cell type and data type

for file in $DIR2/signif_variant*.txt; do
	f=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | cut -d'_' -f3-)
	cmd="$BCFTOOLS view -i "ID=@$file" $VCF -o $DIR_VCF/${f}.vcf.gz -Oz"
	wsbatch --partition=shared-cpu --time=02:00:00 -n1 -c 8 --mem=10000 --wrap="$cmd"
done

 
for file in $DIR_VCF/*.vcf.gz; do
	tabix -p vcf $file
done


#### annotate vcf files

for file in $DIR_VCF/*.vcf.gz; do
	f=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
	$BCFTOOLS annotate -x ^FORMAT/DS -Oz $file -o $DIR_VCF_A/${f}_annotated.vcf.gz
	tabix -p vcf $DIR_VCF_A/${f}_annotated.vcf.gz
	zcat $DIR_VCF_A/${f}_annotated.vcf.gz  | grep -v "##" | cut -f3,10- > $DIR_VCF_A/${f}_annotated2.vcf
done


