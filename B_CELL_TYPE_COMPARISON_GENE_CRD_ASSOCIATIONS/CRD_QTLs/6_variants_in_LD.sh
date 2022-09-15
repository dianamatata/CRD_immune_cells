#!/bin/bash


OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/LD
VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
DIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/conditional_merged
mkdir -p $OUTDIR $OUTDIR/OUT



# compute all the variants of the vcf file
# if LD=1, the variants are in LD, ie in the same haplotype
# here threshold is 0.8, default is 0.9
# are these the variants in LD with r2>0.8 ? what is the threshold exactly, why can't we see the LD value of the pair?

for file in $DIR/*.gz; do
	echo $file
	zcat $file | cut -d' ' -f8 | uniq > $OUTDIR/temp_extract_variants.txt
	cat $OUTDIR/temp_extract_variants.txt | wc -l 
	filename=$( echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
	outfile=$OUTDIR/OUT/${filename}.txt
	cmd="plink --vcf $VCF --extract $OUTDIR/temp_extract_variants.txt --r2 'with-freqs' --ld-window-r2 0.8 --out $OUTDIR/${filename}.txt"
	wsbatch -J plink_$filename \
	--partition=public-cpu \
	--time=00:30:00 -n1 -c 8 \
	--mem=10000 \
	-o OUT/plink_${filename}.out \
	-e OUT/plink_${filename}.err \
	--wrap="$cmd"
done

# temp_extract_variants:         # nbr varies from 32 to more than 70 000


for file in $OUTDIR/*.ld; do
        echo $file
	cat $file | wc -l
done

# constant result: 1300674

<< 'MULTILINE-COMMENT'
        cmd="plink --vcf $VCF --extract $OUTDIR/temp_extract_variants.txt --r2 'with-freqs' --ld-window-r2 0.8 --out $outfile"
        cmd="plink --vcf $VCF --extract $OUTDIR/extract_variants.txt --r2 'with-freqs' --out $OUTDIR/${filename}2.txt"
        cmd="plink --vcf $VCF --extract $OUTDIR/extract_variants.txt --r2 'with-freqs' --ld-window-r2 0.5 --out $OUTDIR/${filename}3.txt"


# compute all the variants linked to this CRD
# maybe better to compute all the variants first
# and then filter for each CRD, otherwise very long
for file in $DIR/*.gz; do
        zcat $file | cut -d' ' -f1 | uniq > $OUTDIR/temp_CRDs.txt
        filename=$( echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
        summary_file=$OUTDIR/${filename}_summary.txt
        echo "CRD variants_in_LD" > $summary_file
        for CRD in $(cat $OUTDIR/temp_CRDs.txt); do
                outfile=$OUTDIR/OUT/${filename}_$CRD.txt
                echo $CRD
                zcat $file | grep $CRD | cut -d' ' -f8 > $OUTDIR/temp_extract_variants.txt
                cmd="plink --vcf $VCF --extract $OUTDIR/temp_extract_variants.txt --r2 'with-freqs' --ld-window-r2 0.8 --out $outfile"
                eval $cmd
                line=echo "$CRD $(echo $outfile | wc -l)"
                cat $line >> $summary_file
        done
done

MULTILINE-COMMENT
