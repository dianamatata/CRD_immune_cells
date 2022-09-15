#### for each case triplet build a nbr_individuals x 3 (triplets) data matrix 
### containing normalized quantifications (individuals in rows, variant, gene and CRD as columns)

# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(qvalue)
library(data.table)
library(tidyverse)
library(dplyr)


#############################################################################################
#
# FUNCTIONS
#
#############################################################################################

keep_shared_CRDs <- function(crd_qtl_cell1_signif,cell1shared){
  # cell1shared=shared_crds$V1
  crd_qtl_cell1_signif_shared=crd_qtl_cell1_signif[(crd_qtl_cell1_signif$phenotype_ID %in% cell1shared), ]
  crd_qtl_cell1_signif_shared
}


#############################################################################################
#
# DIRECTORIES AND DATA
#
#############################################################################################

# directories for mac
dir_vcf='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/vcf'
dir_crd='/Users/dianaavalos/Programming/A_CRD_plots/1_CRD'
dir_rna='/Users/dianaavalos/Programming/Hi-C_correlated_peaks'
dir_rna='/Users/dianaavalos/Programming/A_CRD_plots/RNA'
dir_triplets='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/triplets_signif'
dir_out='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios'
rna_file <- c('/EGAD00001002675_RNA.PC10.bed.gz', '/EGAD00001002674_RNA.PC10.bed.gz', '/EGAD00001002671_RNA.PC10.bed.gz')
names(rna_file) <- c("neut", "mono", "tcell")

# df2 <- as.data.frame(t(vcf[,-1]), stringsAsFactors=F)rectories for cluster
dir_vcf='/home/users/a/avalosma/scratch/12_TRIPLETS/vcf/annotated'
dir_crd='/home/users/a/avalosma/scratch/2_CRD/quantify_ALL'
dir_rna='/home/users/a/avalosma/scratch/4_CRD_residualized'
dir_triplets='/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_signif'
dir_out='/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_data'
rna_file <- c('/EGAD00001002675_RNA.PC10.bed.gz', '/EGAD00001002674_RNA.PC10.bed.gz', '/EGAD00001002671_RNA.PC10.bed.gz')
names(rna_file) <- c("neut", "mono", "tcell")



data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')

for(data_type in data_types){
  for(cell_type in cell_types){
    
    ### loading files
    name_condition=paste0(data_type,'_',cell_type )
    cat(name_condition, '  \n')
    vcf_filtered=fread(paste0(dir_vcf,'/',name_condition,'_annotated2.vcf'), head=TRUE, stringsAsFactors=FALSE)
    vcf=vcf_filtered[grepl("rs",vcf_filtered$ID),] 
    CRD=fread(paste0(dir_crd,'/',name_condition,'.ALLchr.mean.txt.gz'), head=TRUE, stringsAsFactors=FALSE)
    bed=fread(paste0(dir_rna,rna_file[[cell_type]]))
    triplets=fread(paste0(dir_triplets,'/',name_condition,'_triplet.txt'),head=FALSE, stringsAsFactors=FALSE)
    colnames(triplets) = c("variant","gene","CRD")
    triplets$number <- seq.int(nrow(triplets))

        
    ## transpose
    
    vcf2 <- as.data.frame(t(vcf[,-1]), stringsAsFactors=F) 
    vcf2$samples <- colnames(vcf)[-1]
    colnames(vcf2)=vcf$ID # samples in rows, variant in column
    vcf2$samples <- rownames(vcf2)
    
    
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
      t_gene  <- trio$gene
      t_var <- trio$variant
      t_CRD <- trio$CRD
      t_num <- trio$number
      
      if (t_var!='.' && grepl('ENSG', t_gene)){
 
	bed_gene <- bed2[, c(t_gene, "samples")]
	head(bed_gene)
	CRD_CRD <-  CRD2[, c(t_CRD, "samples")]
	head(CRD_CRD)
	vcf_variant <-  vcf2[, c(t_var, "samples")]
	head(vcf_variant)

	merge <- merge(bed_gene,vcf_variant,by.x="samples",by.y="samples",all.x=TRUE)
	merge1 <- merge(merge,CRD_CRD,by.x="samples",by.y="samples",all.x=TRUE)

      	datafr <- na.omit(merge1)
      	triplets_cases[[r]] <- datafr
        write.table(datafr, file = paste0(dir_out,'/',name_condition,'_',t_num,'_',t_var,'_',t_gene,'_',t_CRD, ".txt"), row.names=F, quote=F, col.names=T,sep="\t")

      }
    }
  }
}

# 
### debug
# data_type='hist'
# cell_type='neut'

# ## tryout
# trio <- triplets[1,] 
# t_gene  <- trio$gene
# t_var <- trio$variant
# t_CRD <- trio$CRD
# t_num <- trio$number
# 
# # creation of trio info per sample
# bed_gene <- dplyr::select(bed2,c(t_gene,"samples"))
# head(bed_gene)
# CRD_CRD <- dplyr::select(CRD2,c(t_CRD,"samples"))
# head(CRD_CRD)
# vcf_variant <- dplyr::select(vcf2,c(t_var,"samples"))
# head(vcf_variant)
# 
# merge <- merge(bed_gene,vcf_variant,by.x="samples",by.y="samples",all.x=TRUE)
# merge1 <- merge(merge,CRD_CRD,by.x="samples",by.y="samples",all.x=TRUE)
# 
# datafr <- na.omit(merge1)
# write.table(datafr,file=paste0(dir_out,'/',t_num,'_',t_var,'_',t_gene,'_',t_CRD, ".txt"),row.names=F, quote=F, col.names=T,sep="\t")
# 


