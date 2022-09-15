VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
BCFTOOLS=/srv/beegfs/scratch/groups/funpopgen/Tools/bcftools-1.10.2/bcftools
DIR=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_signif

DIR_VCF=/home/users/a/avalosma/scratch/12_TRIPLETS/vcf_all
DIR_VCF_A=$DIR_VCF/annotated
mkdir -p $DIR_VCF $DIR_VCF_A

#### create vcf files for each file of signif variant per cell type and data type

for file in $DIR/variants_*.txt; do
        f=$(echo $file | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 )
        $BCFTOOLS view -i "ID=@$file" $VCF -o $DIR_VCF/${f}.vcf.gz -Oz
done
