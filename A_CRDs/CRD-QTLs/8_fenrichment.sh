# ---------------------------------------------
# ENRICHMENT ---------------------------------------------
# ---------------------------------------------

INDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/merged 
OUTDIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/significants 

# Step1: Prepare the QTL data  ---------------------------------------------
for file in $OUTDIR/FDR*significant.txt; do
  file=$(echo $file | rev | cut -d "/" -f1 | rev)
  filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '_' -f3-5)
  print $file
  cat $OUTDIR/$file | awk '{ print $9, $10-1, $11, $8, $1, $5 }' | tr " " "\t" | sort -k1,1 -k2,2n > $OUTDIR/$filename.genes.significant.bed
done

# 1. Chromosome ID of the phenotype
# 2. Start position of the phenotype
# 3. End position of the phenotype
# 4. Phenotype ID
# 5. Top variants (not used, can be whatever you want), here CRD
# 6. Strand orientation (important to measure distance between QTLs and phenotypes)


# Step2: Prepare the Phenotype data  ---------------------------------------------
for file in $INDIR/*.txt.gz; do
  file=$(echo $file | rev | cut -d "/" -f1 | rev)
  filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '_' -f1-3)
  print $file
  zcat $INDIR/$file | awk '{ print $2, $3-1, $4, $1, $8, $5 }' | tr " " "\t" | sort -k1,1 -k2,2n > $OUTDIR/$filename.genes.quantified.bed
done



# enrichment PER TFBS  ---------------------------------------------
enrichment_file=/home/users/a/avalosma/scratch/Annotations/TFBS.txt
dir_log=/home/users/a/avalosma/scratch/10_CRD_QTLs/LOG
# only the 50 most present

# /home/users/a/avalosma/scratch/Annotations/TFs_highly_present2.txt
# awk '{print $4}'  TFs_highly_present2.txt | sort | uniq -c # got 50 of them
# awk '{print $4}'  TFs_highly_present2.txt | sort | uniq > list50TFs.txt

dirannot=/home/users/a/avalosma/scratch/Annotations
enrichment=$dirannot/TFBS_sorted_nochr.txt

# for TFBS in $(cat $dirannot/list50TFs.txt)

for TFBS in $(cat $dirannot/all_TFBS_71809153.txt | tail -n +2000) 
do
    echo $TFBS
    cat $enrichment | grep $TFBS > $enrichment_file
    enrichment_file=$dirannot/TFBS_enrichment/${TFBS}_enrichment.txt
    
  for file in $OUTDIR/FDR*significant.txt; do
    file=$(echo $file | rev | cut -d "/" -f1 | rev)
    filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '_' -f3-5)
    cmdTFBS1="QTLtools fenrich --qtl $OUTDIR/$filename.genes.significant.bed \
    --tss $OUTDIR/$filename.genes.quantified.bed \
    --bed $enrichment_file \
    --out $OUTDIR/enrichment_all/${filename}.enrichment.${TFBS}.QTL.in.TF.txt"
  
    wsbatch -J cmdTFBS1.job \
    --partition=public-cpu \
    --time=24:00:00 \
    --mem=10000 \
    -o $dir_log/cmdTFBS1.out \
    -e $dir_log/cmdTFBS1.err \
    --wrap="$cmdTFBS1"
  done
done
    # echo $file





# enrichment epimap  ---------------------------------------------
enrichment_file=/home/users/a/avalosma/scratch/Annotations/Epimap/epimap_2cells.sorted.nochr.bed.gz 
dir_log=/home/users/a/avalosma/scratch/10_CRD_QTLs/LOG

for file in $OUTDIR/FDR*significant.txt; do
  file=$(echo $file | rev | cut -d "/" -f1 | rev)
  filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '_' -f3-5)
  echo $filename
  cmdepimap="QTLtools fenrich --qtl $OUTDIR/$filename.genes.significant.bed \
  --tss $OUTDIR/$filename.genes.quantified.bed \
  --bed $enrichment_file \
  --out $OUTDIR/enrichment/${filename}_epimap.enrichment.QTL.in.TF.txt"
  
  wsbatch -J epimap.job \
  --partition=public-cpu \xs
  --time=24:00:00 \
  --mem=10000 \
  -o $dir_log/epimap.out \
  -e $dir_log/epimap.err \
  --wrap="eval $cmdepimap"
  
done

# enrichment encode  ---------------------------------------------
enrichment_file=/home/users/a/avalosma/scratch/Annotations/encode/encode_3cells_to_22.sorted.nochr.bed.gz
dir_log=/home/users/a/avalosma/scratch/10_CRD_QTLs/LOG


for file in $OUTDIR/FDR*significant.txt; do
  file=$(echo $file | rev | cut -d "/" -f1 | rev)
  filename=$(echo $file | rev | cut -d "/" -f1 | rev | cut -d '_' -f3-5)
  echo $filename
  cmdencode="QTLtools fenrich --qtl $OUTDIR/$filename.genes.significant.bed \
  --tss $OUTDIR/$filename.genes.quantified.bed \
  --bed $enrichment_file \
  --out $OUTDIR/enrichment/${filename}.enrichment.encode.QTL.in.TF.txt"

  wsbatch -J cmdencode.job \
  --partition=public-cpu \
  --time=24:00:00 \
  --mem=10000 \
  -o $dir_log/cmdencode.out \
  -e $dir_log/cmdencode.err \
  --wrap="$cmdencode"
done










# enrichment files concatenate
# need 1 not chr 1
cd ~/scratch/Annotations/encode
# in each subfolder, merge all and sirt and uniq
zcat *.gz >> all.bed
cat all.bed | sort -V -k1,1 -k2,2n | uniq | gzip > encode_monocytes.bed.gz
cat all.bed | sort -V -k1,1 -k2,2n | uniq | gzip > encode_neutrophils.bed.gz
cat all.bed | sort -V -k1,1 -k2,2n | uniq | gzip > encode_tcells.bed.gz

# add column woth origin info
zcat encode_neutrophils.bed.gz | sed "s/$/\tencode_neutrophils/" | gzip > encode_neutrophils1.bed.gz
zcat encode_tcells.bed.gz | sed "s/$/\tencode_tcells/" | gzip > encode_tcells1.bed.gz
zcat encode_monocytes.bed.gz | sed "s/$/\tencode_monocytes/" | gzip > encode_monocytes1.bed.gz

mv encode_neutrophils1.bed.gz encode_neutrophils.bed.gz
mv encode_tcells1.bed.gz encode_tcells.bed.gz
mv encode_monocytes1.bed.gz encode_monocytes.bed.gz

# concatenate
cd /home/users/a/avalosma/scratch/Annotations/encode
zcat neutrophils/encode_neutrophils.bed.gz tcells/encode_tcells.bed.gz monocytes/encode_monocytes.bed.gz >> encode_3cells.bed
cat encode_3cells.bed | sort -V -k1,1 -k2,2n | uniq | gzip > encode_3cells.sorted.bed.gz
zcat /home/users/a/avalosma/scratch/Annotations/encode/encode_3cells.sorted.bed.gz | sed "s/^chr//" | gzip > /home/users/a/avalosma/scratch/Annotations/encode/encode_3cells.sorted.nochr.bed.gz


# epimap   ---------------------------------------------
cd /home/users/a/avalosma/scratch/Annotations/Epimap/CD14_mono
zcat CD14M.bed.gz | sort -V -k1,1 -k2,2n | uniq | gzip > epimap_monocyte.bed.gz
zcat epimap_monocyte.bed.gz | sed "s/$/\tepimap_monocyte/" | gzip > epimap_monocyte1.bed.gz
mv epimap_monocyte1.bed.gz epimap_monocyte.bed.gz

cd /~/scratch/Annotations/Epimap/CD4_Tcells
zcat CD4T.bed.gz | sort -V -k1,1 -k2,2n | uniq | gzip > epimap_tcells.bed.gz
zcat epimap_tcells.bed.gz | sed "s/$/\tepimap_tcells/" | gzip > epimap_tcells1.bed.gz
mv epimap_tcells1.bed.gz epimap_tcells.bed.gz

# concatenate
cd /home/users/a/avalosma/scratch/Annotations/Epimap/CD14_mono
zcat CD4_Tcells/epimap_tcells.bed.gz CD14_mono/epimap_monocyte.bed.gz >> epimap_2cells.bed
cat epimap_2cells.bed | sort -V -k1,1 -k2,2n | uniq | gzip > epimap_2cells.sorted.bed.gz
zcat /home/users/a/avalosma/scratch/Annotations/Epimap/epimap_2cells.sorted.bed.gz | sed "s/^chr//" | gzip \
> /home/users/a/avalosma/scratch/Annotations/Epimap/epimap_2cells.sorted.nochr.bed.gz

# remove weird chr
zcat /home/users/a/avalosma/scratch/Annotations/encode/encode_3cells.sorted.nochr.bed.gz\
|  awk -F "," '$1 <= 22' \
> /home/users/a/avalosma/scratch/Annotations/encode/encode_3cells_filter.bed
cat /home/users/a/avalosma/scratch/Annotations/encode/encode_3cells_filter.bed | gzip > /home/users/a/avalosma/scratch/Annotations/encode/encode_3cells_to_22.sorted.nochr.bed.gz
# JHw09utzbCXRs
# unzip '*.zip'

