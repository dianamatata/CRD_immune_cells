#### for each case triplet build a nbr_individuals x 3 (triplets) data matrix 
### containing normalized quantifications (individuals in rows, variant, gene and CRD as columns)

# Clean environment ---------------------------------
rm(list=ls())
gc()

# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)
library(dplyr)

# Directories ---------------------------------

# Directories for Mac =================================
dir_vcf='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/vcf'
dir_crd='/Users/dianaavalos/Programming/A_CRD_plots/1_CRD'
dir_rna='/Users/dianaavalos/Programming/Hi-C_correlated_peaks'
dir_rna='/Users/dianaavalos/Programming/A_CRD_plots/RNA'
dir_triplets='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/triplets_signif'
dir_out='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/triplets_data'
rna_file <- c('/EGAD00001002675_RNA.PC10.bed.gz', '/EGAD00001002674_RNA.PC10.bed.gz', '/EGAD00001002671_RNA.PC10.bed.gz')
names(rna_file) <- c("neut", "mono", "tcell")

# Directories for Cluster =================================

# df2 <- as.data.frame(t(vcf[,-1]), stringsAsFactors=F)rectories for cluster
dir_vcf='/home/users/a/avalosma/scratch/12_TRIPLETS/vcf/annotated'
dir_crd='/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify_ALL'
dir_rna='/home/users/a/avalosma/scratch/4_CRD_residualized'
dir_triplets='/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_signif'
dir_out='/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_data_shared'
rna_file <- c('/EGAD00001002675_RNA.PC10.bed.gz', '/EGAD00001002674_RNA.PC10.bed.gz', '/EGAD00001002671_RNA.PC10.bed.gz')
names(rna_file) <- c("neut", "mono", "tcell")


# MAIN ---------------------------------

data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
window='1000000'

for(data_type in data_types){
  for(cell_type_ref1 in cell_types){
    for(cell_type_query2 in cell_types){
      if (cell_type_ref1 != cell_type_query2){
        
        ##### first cell is the phenotype, peaks, and 2nd is the CRD location
        name_condition=paste0(data_type,'_',cell_type_query2,'_vs_',cell_type_ref1 )
        cat(name_condition, '  \n')
        
        # Loading Files =================================
        
        vcf_filtered=fread(paste0(dir_vcf,'/',paste0(data_type,'_',cell_type_ref1,'_',window),'_annotated.vcf'), head=TRUE, stringsAsFactors=FALSE)
        vcf=vcf_filtered[grepl("rs",vcf_filtered$ID),] 
        
        # CRD phenotype in cell2 and gene exp in cell 2
        CRD=fread(paste0(dir_crd,'/',name_condition,'.ALLchr.mean.txt.gz'), head=TRUE, stringsAsFactors=FALSE)
        bed=fread(paste0(dir_rna,rna_file[[cell_type_query2]]))
        
        triplets=fread(paste0(dir_triplets,'/',paste0(data_type,'_',cell_type_ref1,'_',window),'_triplet.txt'),head=FALSE, stringsAsFactors=FALSE)

        colnames(triplets) = c("variant","gene","CRD", "dist")
        triplets$number <- seq.int(nrow(triplets))
        
        
        ## transpose
        
        vcf2 <- as.data.frame(t(vcf[,-(1:9)]), stringsAsFactors=F) 
        vcf2$samples <- colnames(vcf[,-(1:9)])
        colnames(vcf2)=vcf$ID # samples in rows, variant in column
        vcf2$samples <- rownames(vcf2)
        vcf$INFO <- as.double(gsub("AF=", "", vcf$INFO))   
        
        
        bed2 <- as.data.frame(t(bed[,-c(1:6)]), stringsAsFactors=F) 
        colnames(bed2)=bed$id # gene_list in col
        rownames(bed2)=colnames(bed[,-c(1:6)]) # samples in row
        bed2$samples <- rownames(bed2)
        
        
        rownames(CRD) <- CRD$id
        CRD2 <- as.data.frame(t(CRD[,-c(1:6)]), stringsAsFactors=F) 
        colnames(CRD2)=CRD$id # crd id in col
        CRD2$samples <- rownames(CRD2)
        
        ## loop over all triplets
        triplets_cases <- list()
        for(r in 1:nrow(triplets)){
          trio <- triplets[r,] 
          t_var <- trio$variant
          t_gene  <- trio$gene
          t_CRD <- trio$CRD
          t_num <- trio$number
          t_dist <- trio$dist
          
          # cat(t_var, ' ',t_gene,' ',t_CRD,' ',t_num )
          
          if (!(t_gene %in% colnames(bed2)) ){
            write(paste0(t_gene,'  '),file=paste0(dir_out,'/',name_condition,'_genes_out.txt',append=TRUE))
          }
          
          if (t_var!='.' && grepl('ENSG', t_gene) && t_gene %in% colnames(bed2)){
            
            bed_gene <- bed2[, c(t_gene, "samples")]
            head(bed_gene)
            
            CRD_CRD <-  CRD2[, c(t_CRD, "samples")]
            head(CRD_CRD)
            
            vcf_variant <-  vcf2[, c(t_var, "samples")]
            vcf_variant[,1]=as.double(substr(vcf_variant[,1],1,1))
            head(vcf_variant)
            
            
            variant_AF= vcf %>% filter(ID==t_var)
            variant_AF=variant_AF$INFO
            
            merge <- merge(bed_gene,vcf_variant,by.x="samples",by.y="samples",all.x=TRUE)
            merge1 <- merge(merge,CRD_CRD,by.x="samples",by.y="samples",all.x=TRUE)
            
            datafr <- na.omit(merge1)
            triplets_cases[[r]] <- datafr
            write.table(datafr, file = paste0(dir_out,'/',name_condition,'_',t_num,'_',t_var,'_',t_gene,'_',t_CRD,'_',t_dist, '_', toString(variant_AF), ".txt"), row.names=F, quote=F, col.names=T,sep="\t")
          } 
        }
      }
    }
  }
}


# 
### debug
data_type='hist'
cell_type_ref1='neut'
cell_type_query2='mono'

