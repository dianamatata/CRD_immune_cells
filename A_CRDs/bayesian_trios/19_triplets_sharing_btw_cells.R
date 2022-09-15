
# Clean environment ---------------------------------
rm(list=ls())
gc()


# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)



# FUNCTION ---------------------------------

add_colnames_signif <- function(triplets_cell1_signif){
  colnames(triplets_cell1_signif) = c("phenotype_ID","CRD","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                      "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end","rank",
                                      "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
  triplets_cell1_signif
}


add_colnames_all <- function(triplets_cell2_all){
  colnames(triplets_cell2_all) = c("phenotype_ID","CRD","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                   "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end",
                                   "nom_pval","r_squared","slope","besthit")
  triplets_cell2_all
}


keep_shared_CRDs <- function(triplets_cell1_signif,cell1shared){
  triplets_cell1_signif_shared=triplets_cell1_signif[(triplets_cell1_signif$CRD %in% cell1shared), ]
  triplets_cell1_signif_shared
}

compute_triplets_sharing <- function(shared_crds,triplets_cell1_signif_shared,triplets_cell2_all_shared,out_directory, name){
  
  filename=paste0(out_directory,name,'_pvals_if_QTL_CRD_shared.txt')
  file.create(filename)
  
  for (i in 1:dim(triplets_cell1_signif_shared)[1]) {
    # cat(i,' ')
    
    # find equiv CRDs in cell2 + variant
    cell2_equivalent_crd1 = shared_crds$V2[which(shared_crds$V1 == triplets_cell1_signif_shared[i,]$CRD)] # many CRDs can correspond
    # filter
    subset=triplets_cell2_all_shared %>% filter(CRD %in%  cell2_equivalent_crd1) %>% 
      filter(var_ID %in% triplets_cell1_signif_shared[i,]$var_ID) %>% 
      filter(phenotype_ID %in% triplets_cell1_signif_shared[i,]$phenotype_ID)
    
    subset$nom_pval
    
    # save values
    
    if (length(subset[,1]) != 0){
      line=paste0(subset$nom_pval)
      line_bis=paste0(subset$nom_pval, ' ',subset$var_ID, ' ',subset$CRD)
      write(line_bis,file=filename,append=TRUE)
    }
  }
  
  result=as.data.frame(data.table::fread(filename), head=FALSE, stringsAsFactors=FALSE)
  result 
}


# Directories ---------------------------------

# for cluster

sharing_directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/CRD_sharing/'
# and filter nominal with significant CRDs, and genes?
dir_triplets_all='/home/users/a/avalosma/scratch/12_TRIPLETS/not_signif/merged_nominal1/triplets/'
dir_triplets_signif='/home/users/a/avalosma/scratch/12_TRIPLETS/significants/triplets/'
out_directory='/home/users/a/avalosma/scratch/12_TRIPLETS/'

args <- commandArgs(trailingOnly = TRUE)
data_type <- args[1]
module='mean'
cell1 <- args[2]
cell2 <- args[3]


# MAIN ---------------------------------

FDRthreshold=0.05
cell_pairs=c(c('neut','mono'),c('neut','tcell'),c('mono','neut'),c('mono','tcell'),c('tcell','neut'),c('tcell','mono'))

name=paste0(data_type,'_',module,'_',cell1,'_vs_',cell2)
cat (name, '  ')

triplets_cell1_signif=fread(paste0(dir_triplets_signif,'FDR_',FDRthreshold,'_',data_type,'_',cell1,'_CRD_gene_var_ALL.significant.txt'), head=FALSE,stringsAsFactors=FALSE,fill=TRUE)
triplets_cell2_all=fread(paste0(dir_triplets_all,data_type,'_',cell2,'_CRD_gene_var_ALL.txt'), head=FALSE, stringsAsFactors=FALSE)

shared_crds=as.data.frame(data.table::fread(paste0(sharing_directory,name,'_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE))

triplets_cell1_signif=add_colnames_signif(triplets_cell1_signif)
triplets_cell2_all=add_colnames_all(triplets_cell2_all) 

triplets_cell1_signif_shared=keep_shared_CRDs(triplets_cell1_signif,shared_crds$V1)
triplets_cell2_all_shared=keep_shared_CRDs(triplets_cell2_all,shared_crds$V2)

# check sizes
cat( ' triplets_cell1_signif_shared ', dim(triplets_cell1_signif_shared)[1] ,
     ' triplets_cell1_signif ', dim(triplets_cell1_signif)[1],
     ' triplets_cell2_all_shared ', dim(triplets_cell2_all_shared)[1] ,
     ' triplets_cell2_all ', dim(triplets_cell2_all)[1] ,
     ' shared_crds ', dim(shared_crds)[1] )

hist(triplets_cell2_all_shared$nom_pval)
length(unique(triplets_cell2_all_shared$var_ID))
# compute sharing
res=compute_triplets_sharing(shared_crds,triplets_cell1_signif_shared,triplets_cell2_all_shared,out_directory, name)

# plots
qobj=qvalue(res$V1)
pi0=qobj$pi0
pi1=1-pi0
pdf(paste0(plot_directory,"histogram_p_",name,".pdf"))
hist(res$V3, main=paste0(' pi0=',round(pi0,digits = 3),'  pi1=',round(pi1,digits = 3)), xlab='p values', breaks=20,cex.main=1.5,  cex.lab=1.2,cex.axis=1.2)
dev.off()


