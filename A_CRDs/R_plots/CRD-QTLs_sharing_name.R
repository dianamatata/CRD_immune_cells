### goal: compute pi1 estimate of CRD-QTL sharing

# we need the merged and the 5%FDR significant for each
# shared file

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(qvalue)
library(ggplot2)
library(tidyverse)
library(qvalue)

#############################################################################################
#
# FUNCTION
#
#############################################################################################

add_colnames_signif <- function(crd_qtl_cell1_signif){
  
  colnames(crd_qtl_cell1_signif) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                   "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end","rank",
                                   "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
crd_qtl_cell1_signif
}

add_colnames_all <- function(crd_qtl_cell2_all){
  colnames(crd_qtl_cell2_all) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                         "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end",
                         "nom_pval","r_squared","slope","besthit")
  
  # colnames(crd_qtl_cell2_all) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
  #                                 "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end","dof1",
  #                                 "dof2","bml1","bml2","nom_pval","r_squared","slope","slope_se","adj_pvalue")
  crd_qtl_cell2_all
}


keep_shared_CRDs <- function(crd_qtl_cell1_signif,cell1shared){
  # cell1shared=shared_crds$V1
  crd_qtl_cell1_signif_shared=crd_qtl_cell1_signif[(crd_qtl_cell1_signif$phenotype_ID %in% cell1shared), ]
  crd_qtl_cell1_signif_shared
}

compute_CRD_QTL_sharing <- function(shared_crds,crd_qtl_cell1_signif_shared,crd_qtl_cell2_all_shared,out_directory, name){
  
  filename=paste0(out_directory,name,'_pvals_if_QTL_CRD_shared.txt')
  file.create(filename)
  
  for (i in 1:length(crd_qtl_cell1_signif_shared[,1])) {
    # cat(i,' ')

    # find equiv CRDs in cell2 + variant
    cell2_equivalent_crd1 = shared_crds$V2[which(shared_crds$V1 == crd_qtl_cell1_signif_shared[i,]$phenotype_ID)] # many CRDs can correspond
    variant = crd_qtl_cell1_signif_shared[i,]$var_ID
    # filter
    sub1=crd_qtl_cell2_all_shared %>% filter(phenotype_ID %in%  cell2_equivalent_crd1)   
    subset=sub1 %>% filter(var_ID %in%  variant) 
    subset$nom_pval
    
    # save values
    
    if (length(subset[,1]) != 0){
      line=paste0(subset$nom_pval)
      line_bis=paste0(subset$nom_pval, ' ',subset$var_ID, ' ',subset$phenotype_ID)
      write(line_bis,file=filename,append=TRUE)
    }
  }
  
  result=as.data.frame(data.table::fread(filename), head=FALSE, stringsAsFactors=FALSE)
  result 
}



#############################################################################################
#
# MAIN
#
#############################################################################################
#### for cluster

directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/CRD_sharing/'
out_directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/shared/'
dir_CRD_QTL_all='/home/users/a/avalosma/scratch/10_CRD_QTLs/merged_nominal_no_FDR/'
dir_CRD_QTL_signif='/home/users/a/avalosma/scratch/10_CRD_QTLs/conditional_merged/'
# file_est=paste0(out_directory,name,'_pi_estimates.txt')
# file.create(file_est)

cell_pairs=c(c('neut','mono'),c('neut','tcell'),c('mono','neut'),c('mono','tcell'),c('tcell','neut'),c('tcell','mono'))
FDRthreshold=0.05
module='mean'
data_type='methyl'
i=1

args <- commandArgs(trailingOnly = TRUE)
data_type <- args[1]
module <- args[2]
cell1 <- args[3]
cell2 <- args[4]


name=paste0(data_type,'_',module,'_',cell1,'_vs_',cell2)
cat (name, '  ')

#crd_qtl_cell1_signif=as.data.frame(data.table::fread(paste0(dir_CRD_QTL_signif,'FDR_',FDRthreshold,'_',data_type,'_',cell1,'_',module, '_ALL.significant.txt'), head=FALSE, stringsAsFactors=FALSE))
crd_qtl_cell1_signif=as.data.frame(data.table::fread(paste0(dir_CRD_QTL_signif,data_type,'_',cell1,'_',module,'__ALL.txt.gz'), head=FALSE, stringsAsFactors=FALSE))
crd_qtl_cell2_all=as.data.frame(data.table::fread(paste0(dir_CRD_QTL_all,data_type,'_',cell2,'_',module,'__ALL.txt.gz'), head=FALSE, stringsAsFactors=FALSE))
shared_crds=as.data.frame(data.table::fread(paste0(directory,name,'_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE))

crd_qtl_cell1_signif=add_colnames_signif(crd_qtl_cell1_signif)
crd_qtl_cell2_all=add_colnames_all(crd_qtl_cell2_all)

crd_qtl_cell1_signif_shared=keep_shared_CRDs(crd_qtl_cell1_signif,shared_crds$V1)
crd_qtl_cell2_all_shared=keep_shared_CRDs(crd_qtl_cell2_all,shared_crds$V2)

# check sizes
cat( ' crd_qtl_cell1_signif_shared ', length(crd_qtl_cell1_signif_shared[,1]) ,
     ' crd_qtl_cell1_signif ', length(crd_qtl_cell1_signif[,1]),
     ' crd_qtl_cell2_all_shared ', length(crd_qtl_cell2_all_shared[,1]),
     ' crd_qtl_cell2_all ', length(crd_qtl_cell2_all[,1]),
     ' shared_crds ', length(shared_crds[,1]))

hist(crd_qtl_cell2_all_shared$nom_pval)
length(unique(crd_qtl_cell2_all_shared$var_ID))
# compute sharing
res=compute_CRD_QTL_sharing(shared_crds,crd_qtl_cell1_signif_shared,crd_qtl_cell2_all_shared,out_directory, name)

# plots
qobj=qvalue(res$V1)
pi0=qobj$pi0
pi1=1-pi0
pdf(paste0(plot_directory,"histogram_p_",name,".pdf"))
hist(res$V3, main=paste0(' pi0=',round(pi0,digits = 3),'  pi1=',round(pi1,digits = 3)), xlab='p values', breaks=20,cex.main=1.5,  cex.lab=1.2,cex.axis=1.2)
dev.off()





