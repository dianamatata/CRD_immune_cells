# CRD paper
# Genetic Variants Associated with Blood Count Phenotypes

# -------------------------------------------
# all the Chen QTLs, FDR is column 7 #extract signif at FDR 5%

DIR_Chen_QTLs=/home/users/a/avalosma/scratch/external_datasets/Chen2016/qtl_as/QTL_RESULTS
DIR_Chen_QTLs_signif=/home/users/a/avalosma/scratch/external_datasets/Chen2016/qtl_as/QTL_RESULTS_signif
DIR_CRD_QTLcis=/home/users/a/avalosma/scratch/CRD_project_outputs/10_CRD_QTLs/significants
DIR_CRD_QTLtrans_aCRD=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0/TRANS-EQTL/Q0.01/aCRD
DIR_CRD_QTLtrans_sCRD=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0/TRANS-EQTL/Q0.01/sCRD
DIR_QTL=/home/users/a/avalosma/scratch/CRD_project_outputs/QTL_list
DIR_LD=/home/users/a/avalosma/scratch/CRD_project_outputs/LD_QTL_list
# CADD
DIR_CADD=/home/users/a/avalosma/scratch/external_datasets/CADD_37 
file_CADD=whole_genome_SNVs.tsv.gz
# VCF for LD
VCF_all_Blueprint=/home/users/a/avalosma/scratch/CRD_project_outputs/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz

# GWAS_files diseases
DIR_GWAS=/home/users/a/avalosma/scratch/external_datasets/GWAS_diseases
GWAS_UC=$DIR_GWAS/UC/26192919-GCST003045-EFO_0000729-Build37.f.tsv.gz
GWAS_DT1=$DIR_GWAS/DT1/25751624-GCST005536-EFO_0001359-Build37.f.tsv.gz
GWAS_RA=$DIR_GWAS/RA/24390342-GCST002318-EFO_0000685-Build37.f.tsv.gz
GWAS_CEL=$DIR_GWAS/CEL/20190752-GCST000612-EFO_0001060-Build37.f.tsv.gz
GWAS_CD=$DIR_GWAS/CD/26192919-GCST003044-EFO_0000384-Build37.f.tsv.gz
GWAS_IBD=$DIR_GWAS/IBD/26192919-GCST003043-EFO_0003767-Build37.f.tsv.gz
GWAS_MS=$DIR_GWAS/MS/24076602-GCST005531-EFO_0003885-Build37.f.tsv.gz
GWAS_DT2_36=$DIR_GWAS/DT2/22885922-GCST005047-EFO_0001360-Build36.f.tsv.gz

# DIR_OUT=/home/users/a/avalosma/scratch/CRD_project_outputs/CADD_QTLs


# get signif QTLs chen
for f in $DIR_Chen_QTLs/*
do
	base_name=$(basename ${f})
	echo "Processing $base_name"
	zcat $f | awk '{ if($7 <= 0.05) { print }}' | tr ' ' '\t'  | gzip > $DIR_Chen_QTLs_signif/$base_name 
done
# weird because lots of non signif # ASK DIOGO

# Get all QTLs Chen: have the chr, pos, ref and alt alleles separated for CADD
for f in $DIR_Chen_QTLs_signif/*all_summary.txt.gz
do
	name=$(echo $(basename ${f}) | sed 's/_all_summary.txt.gz//')
	echo "Processing $name"
	zcat $f | tr ':' '\t' | tr '_' '\t' | bgzip > tempfile
	mv tempfile $f
done

# -------------------------------------------
# Get ALL significant QTLs, only RS values

 # create list of Chen QTLs to test with LD
for f in $DIR_Chen_QTLs_signif/*all_summary.txt.gz
do
	name=$(echo $(basename ${f}) | sed 's/_all_summary.txt.gz//')
	echo "Processing $name"
	zcat $f  | cut -f5 > $DIR_QTL/SNPlist_Chen_${name}.txt
done

# Get all CRDs-QTL cis
for f in $DIR_CRD_QTLcis/*significant.txt
do
	name=$(echo $(basename ${f}) | sed 's/FDR_0.05_//' | sed 's/__ALL.significant.txt//')
	echo "Processing $name"
	cat $f | cut -f8 -d' ' > $DIR_QTL/SNPlist_cisCRDQTLs_${name}.txt
done

# Get all CRDs-QTL trans
for f in $DIR_CRD_QTLtrans_aCRD/*significant.txt
do
	name=$(echo $(basename ${f}) | sed 's/.significant.txt//' | sed 's/transeqtl_//')
	echo "Processing $name"
	cat $f | cut -f6 -d' ' > $DIR_QTL/SNPlist_trans_aCRDQTLs_${name}.txt
done

for f in $DIR_CRD_QTLtrans_sCRD/*significant.txt
do
	name=$(echo $(basename ${f}) | sed 's/.significant.txt//' | sed 's/transeqtl_//')
	echo "Processing $name"
	cat $f | cut -f6 -d' ' > $DIR_QTL/SNPlist_trans_sCRDQTLs_${name}.txt
done

## 1) get SNPs that are in LD in 500kb window r2 >= 0.8 for variants in QTLs

for f in $DIR_QTL/*SNPlist_*
do
	name=$(echo $(basename ${f}) | sed 's/.txt//' | sed 's/SNPlist_//')
	echo $name
	if [ -s $f ]
	then
     echo "File not empty"
     plink --vcf $VCF_all_Blueprint \
	  --extract $f \
	  --r2 'with-freqs' \
	  --ld-window-kb 500  \
	  --ld-window-r2 0.8 \
	  --out $DIR_LD/${name}.txt
	else
     echo "File empty"
	fi
done

## 2) Overlap with GWAS: extract all the SNPs in our list that overlap with GWAS hits
# get all GWAS ordered
for GWAS_file in $GWAS_UC $GWAS_CD $GWAS_IBD $GWAS_RA $GWAS_CEL $GWAS_MS $GWAS_DT1 $GWAS_DT2_36
do
	name=$(echo $(basename ${GWAS_file}) | sed 's/-Build37.f.tsv.gz//')
	echo $(basename ${name})
	zcat $GWAS_UC | awk '{print($3,"\t",$5,"\t",$1,"\t",$2)}' > $DIR_overlap/GWAS_${name}.txt
done

# have all variants in list
for f in $DIR_LD/trans*.txt.ld
do
	name=$(echo $(basename ${f}))
	echo $name
	cat $f | sed 's/  */\t/g' | cut -f2,3 > $DIR_overlap/${name}.txt
	cat $f | sed 's/  */\t/g' | cut -f6,7 >> $DIR_overlap/${name}.txt
done

# TODO:  then in R compute overlap, easier than bash.....
# get pvalues

# GWAS blood traits
GWAS_sum_stats_bloodtraits=/home/users/a/avalosma/scratch/external_datasets/Vuckovic2020/UKBB_blood_cell_traits
for blood_trait in plt.assoc rbc.assoc ret.assoc mono.assoc neut.assoc eo.assoc baso.assoc lymph.assoc wbc.assoc
do
	echo $blood_trait
	file=$GWAS_sum_stats_bloodtraits/$blood_trait
	cat $file | sed 's/,/\t/g' |  awk '{print($3,"\t",$5,"\t",$2,"\t",$10)}'
done
# TODO:  then in R compute overlap or gWAS blood, easier than bash.....




# -------------------------------------------
# Apply CADD score to list of variants or vcf
# apply CADD score to vcf!!

# filter all variants of QTL, 3	135120434	rs1069507	C	T
bcftools query -f '%CHROM %POS %REF %ALT\n' $VCF_all_Blueprint > $DIR_CADD/VCF_all_Blueprint_4cols.txt
# 1) apply CADD score to $DIR_CADD/VCF_all_Blueprint_4cols.txt

# 2) in $DIR_CADD/VCF_all_Blueprint_4cols.txt extract all rs values of qtlfile, extract PHRED column

zcat $VCF_all_Blueprint | grep 'rs1069507'
# 1	14376994	rs35640915	T	G
# 1	14377073	rs726596	G	A

cat $DIR_CADD/CADD_sample.txt | head
#Chrom	Pos	Ref	Alt	RawScore	PHRED
# 1	10001	T	A	0.088260	4.066
# 1	10001	T	C	0.105475	4.













/srv/beegfs/scratch/users/a/avalosma/external_datasets/CADD_37
/Users/dianaavalos/Desktop/reviews avalos/GWAS_overlap/

https://github.com/kircherlab/CADD-scripts


/Users/dianaavalos/Desktop/reviews\ avalos/CRD_immune_cells-main/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/DIR_CRD_QTLs/6_variants_in_LD.sh















# compute overlap
# GARFIELD
DIR_GARFIELD=/home/users/a/avalosma/bin/GARFIELD
DIR_G_OUT=/home/users/a/avalosma/scratch/CRD_project_outputs/GARFIELD_output
chr=22
prunetags=$DIR_GARFIELD/garfield-data/tags/r01/chr${chr}
clumptags=$DIR_GARFIELD/garfield-data/tags/r08/chr${chr}
maftssd=$DIR_GARFIELD/garfield-data/maftssd/chr${chr}
pvalue=$DIR_GARFIELD/garfield-data/pval/GIANT_HEIGHT/chr${chr}
annot=$DIR_GARFIELD/garfield-data/annotation/chr${chr}
linkfile=$DIR_GARFIELD/garfield-data/annotation/link_file.txt
pthreshs=1e-5
binning=n5,m5,t5
# -c specifies if GARFIELD is to run one annotation at a time (0) or if additionally heuristic model selection is to be performed (1). 
#  Note: to use -c 1 GARFIELD must first be run with the -c 0 option
# Single annot level enrichment (0) or conditional analysis (1) # Default is 0, takes all annots in link file
thresh=0.05
min_variants_per_tresh=10

prepfile=$DIR_G_OUT/test1
outfile=$DIR_G_OUT/test1.out
outfile2=$DIR_G_OUT/test1.out2
outfile2_1=$DIR_G_OUT/test1.out2_c1

output_path_prefix=$DIR_G_OUT/test1_
plot_title=test1

# prune_and_clump
# prepares the input data for garfield-test.R and garfield-Meff-Padj.R.
$DIR_GARFIELD/garfield-v2/garfield-prep-chr -ptags $prunetags -ctags $clumptags -maftss $maftssd -pval $pvalue -ann $annot -o $prepfile -chr $chr

# computes the effective number of annotations and adjusted for multiple testing p-values from the output of the garfield-prep-chr tool.
Rscript $DIR_GARFIELD/garfield-v2/garfield-Meff-Padj.R -i $prepfile -o $outfile

# computes odds ratios and enrichment p-values, while using accounting for variant differences in number of high LD proxies (clump_tag_count), 
# minor allele frequency (MAF) and transcription start site (TSS) distance. Variants for each of these three features are split into custom number of bins
# and are represented by categorical covariates in a logistic regression model. condition=0
Rscript $DIR_GARFIELD/garfield-v2/garfield-test.R  -i $prepfile -o $outfile2 -l $linkfile -pt 1,0.1,0.01,1e-4,1e-5,1e-6,1e-7,1e-8 -b $binning \
 -c 

Rscript $DIR_GARFIELD/garfield-v2/garfield-test.R  -i $prepfile -o $outfile2_1 -l $linkfile -pt 1,0.1,0.01,1e-4,1e-5,1e-6,1e-7,1e-8 -b $binning \
 -c 1 -padj 0.01
 # [-s subset] [-ct condthresh -padj padj]
 # cat test1.out2 | cut -f2 -d' ' | less shows all thresholds used

# produces the final figures from the table of results. 
Rscript $DIR_GARFIELD/garfield-v2/garfield-plot.R -i $outfile2 -o $output_path_prefix -t $plot_title -f $min_variants_per_tresh  -padj $thresh 
-c "black,grey,red,pink,blue,white,green,yellow"

# deepskyblue4, deepskyblue, aquamarine3, coral, darkorange, bisque2”
 # -s subset [-col set_of_colour_names]
 # which values???
Rscript $DIR_GARFIELD/garfield-v2/garfield-plot.R -i $outfile2 -o $DIR_G_OUT/f10 -t f10 -f 10  -padj  0.05
Rscript $DIR_GARFIELD/garfield-v2/garfield-plot.R -i $outfile2 -o $DIR_G_OUT/f10_0.05 -t f10 -f 10  -padj  0.005




# interesct with R
# m is in the matched string # https://stackoverflow.com/questions/67314459/intersecting-to-files-by-several-columns-using-awk
cmd="awk 'BEGIN{FS=OFS=","} NR==1{print "Chrom\tPos\tRef\tAlt\tPHRED"} FNR==NR {m1[$1]=$1; m1[$2]=$2; m1[$3]=$3; m1[$4]=$4;next} FNR>1 {print $0, m1[$1], m1[$2], m2[$1], m2[$2]}' $test_file $DIR_CADD/$file_CADD > $DIR_QTL/test_intersect.txt"

# otherwise test_original_plink_cohorts_common_SNPs.R





OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/LD
VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/
DIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/conditional_merged
mkdir -p $OUTDIR $OUTDIR/OUT



# 29 blood cell phenotypes
GWAS_sum_stats_bloodtraits=/home/users/a/avalosma/scratch/CRD_project_outputs/Vuckovic2020 

# trait-variant associations
# normal correlation? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4423381/



awk 'BEGIN{FS=OFS=","} NR==1{print ""} FNR==NR {m1[$1]=$1; m1[$2]=$2;next} FNR>1 {print $0, m1[$1], m1[$2]}' $DIR_overlap/templist.txt  $DIR_overlap/GWAS_all.txt  > $DIR_overlap/test_intersect.txt

awk 'BEGIN{FS=OFS=","} NR==1{print ""} FNR==NR {m1[$1]=$1; m1[$2]=$2;next} FNR>1 {print m1[$1], m1[$2]}'  $DIR_overlap/GWAS_all.txt $DIR_overlap/templist.txt  > $DIR_overlap/test_intersect.txt
cat $DIR_overlap/test_intersect.txt | head







