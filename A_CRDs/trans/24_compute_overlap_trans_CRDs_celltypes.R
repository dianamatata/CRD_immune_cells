# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(tidyverse)
library(data.frame)
library(R.utils) # count lines
library(corrplot)



#############################################################################################
#
# FUNCTION
#
#############################################################################################



compute_shared_transCRD_v2 <- function(shared_crds,trans_crd_cell1_signif,trans_crd_cell2_all,filename){
  # function written for cell1_vs_cell2_sharedCRDs.txt, applicable for all pairs afterwards
  
  # have ALL with only the shared part, save the size of these
  
  # have significant with only the shared part
  
  # for the significant cell1, look for the p value in all cell 2
  
  # compute the pi1 estimate
  
  
  # select CRD IDs in significant trans associations
  trans_crd_cell1_CRD_pairs=trans_crd_cell1 %>% select(4,8)
  trans_crd_cell2_CRD_pairs=trans_crd_cell2 %>% select(4,8)
  
  # select PAIRS in trans_crd_cell1_CRD_pairs that are also in shared_crds
  cell1shared=shared_crds$V1
  cell2shared=shared_crds$V2
  
  trans_crd_cell1_CRD_pairs_shared=trans_crd_cell1_CRD_pairs[(trans_crd_cell1_CRD_pairs$V4 %in% cell1shared) & (trans_crd_cell1_CRD_pairs$V8 %in% cell1shared), ]
  trans_crd_cell2_CRD_pairs_shared=trans_crd_cell2_CRD_pairs[(trans_crd_cell2_CRD_pairs$V4 %in% cell2shared) & (trans_crd_cell2_CRD_pairs$V8 %in% cell2shared), ]
  
  file.create(filename)
  # replace in trans_crd_cell1_CRD_pairs_shared with the equivalent in cell2 -B
  for (i in 1:length(trans_crd_cell1_CRD_pairs_shared[,1])) {
    crd_pair_cell1=c(trans_crd_cell1_CRD_pairs_shared[i,]$V4,trans_crd_cell1_CRD_pairs_shared[i,]$V8)
    cell2_equivalent_crd1 = shared_crds$V2[which(shared_crds$V1 == trans_crd_cell1_CRD_pairs_shared[i,]$V4)]
    cell2_equivalent_crd2 = shared_crds$V2[which(shared_crds$V1 == trans_crd_cell1_CRD_pairs_shared[i,]$V8)]
    crd_pair_cell1_in_cell2=c(cell2_equivalent_crd1 ,cell2_equivalent_crd2) # equivalent of cell1_pair in cell2_pair
    # pair_test=c("10_internal_7905",  "19_internal_6283") # idx  5160 : 10_internal_7905  19_internal_6283, from trans_crd_cell2_CRD_pairs_shared
    # find cell2_pair in trans_crd_cell2_CRD_pairs_shared
    indexes1=which(trans_crd_cell2_CRD_pairs_shared$V8 == crd_pair_cell1_in_cell2[1] &  trans_crd_cell2_CRD_pairs_shared$V4 == crd_pair_cell1_in_cell2[2])
    indexes2=which(trans_crd_cell2_CRD_pairs_shared$V4 == crd_pair_cell1_in_cell2[1] &  trans_crd_cell2_CRD_pairs_shared$V8 == crd_pair_cell1_in_cell2[2])
    for (index in c(indexes1,indexes2)){
      if (length(trans_crd_cell2_CRD_pairs_shared[index,1]) != 0){
        line=paste0(trans_crd_cell2_CRD_pairs_shared[index,]$V4,' ',trans_crd_cell2_CRD_pairs_shared[index,]$V8)
        write(line,file=filename,append=TRUE)
      }
    }
  }
  result=as.data.frame(data.table::fread(filename), head=FALSE, stringsAsFactors=FALSE)
  result 
}

#############################################################################################
compute_CRD_QTL_sharing <- function(shared_crds,trans_crd_cell1_signif,trans_crd_cell2_all,filename){

  trans_crd_cell1_signif_shared <- data.frame(col_name=character(0),unique_cnt=integer(0))
  
  # select PAIRS that are also in shared_crds
  trans_crd_signif_cell1_shared=trans_crd_cell1_signif[(trans_crd_cell1_signif$id1 %in% shared_crds$V1) & (trans_crd_cell1_signif$id2 %in% shared_crds$V1), ]
  trans_crd_all_cell2_shared=trans_crd_cell2_all[(trans_crd_cell2_all$id1 %in% shared_crds$V2) & (trans_crd_cell2_all$id2 %in% shared_crds$V2), ]
  
  # add id1_equiv, id2_equiv
  for (i in 1:length(trans_crd_signif_cell1_shared[,1])) {
    # sometimes many results
    trans_crd_signif_cell1_shared[i,]$id1_equiv=shared_crds$V2[which(shared_crds$V1 == trans_crd_signif_cell1_shared[i,]$id1)]
  }
  
  cat( ' trans_crd_cell1_signif ', length(trans_crd_cell1_signif[,1]),
       ' trans_crd_signif_cell1_shared ', length(trans_crd_signif_cell1_shared[,1]),
       ' trans_crd_cell2_all ', length(trans_crd_cell2_all[,1]),
       ' trans_crd_all_cell2_shared ', length(trans_crd_all_cell2_shared[,1]),
       ' shared_crds ', length(shared_crds[,1]))
  

  # find cell1_pair in trans_crd_cell2_CRD_pairs_shared
  for (i in 1:length(trans_crd_signif_cell1_shared[,1])) {
    
    crd_pair_cell1_in_cell2=c(trans_crd_signif_cell1_shared$id1,trans_crd_signif_cell1_shared$id2)
    
    indexes1=which(trans_crd_all_cell2_shared$id2 == crd_pair_cell1_in_cell2[1] &  trans_crd_all_cell2_shared$id1 == crd_pair_cell1_in_cell2[2])
    indexes2=which(trans_crd_all_cell2_shared$id1 == crd_pair_cell1_in_cell2[1] &  trans_crd_all_cell2_shared$id2 == crd_pair_cell1_in_cell2[2])
    
    for (index in c(indexes1,indexes2)){
      if (length(trans_crd_cell2_CRD_pairs_shared[index,1]) != 0){
        line=paste0(trans_crd_cell2_CRD_pairs_shared[index,]$adj_pval)
        row=trans_crd_cell2_CRD_pairs_shared[index,]
        crd_qtl_cell1_signif_shared=rbind(crd_qtl_cell1_signif_shared,row)
        write(line,file=filename,append=TRUE)
      }
    }
  }
  
}
  

#############################################################################################



compute_shared_transCRD <- function(shared_crds,trans_crd_cell1,trans_crd_cell2,filename){
  # function written for cell1_vs_cell2_sharedCRDs.txt, applicable for all pairs afterwards
  
  
  # select CRD IDs in significant trans associations
  trans_crd_cell1_CRD_pairs=trans_crd_cell1 %>% select(4,8)
  trans_crd_cell2_CRD_pairs=trans_crd_cell2 %>% select(4,8)
  
  # select PAIRS in trans_crd_cell1_CRD_pairs that are also in shared_crds
  cell1shared=shared_crds$V1
  cell2shared=shared_crds$V2
  
  trans_crd_cell1_CRD_pairs_shared=trans_crd_cell1_CRD_pairs[(trans_crd_cell1_CRD_pairs$V4 %in% cell1shared) & (trans_crd_cell1_CRD_pairs$V8 %in% cell1shared), ]
  trans_crd_cell2_CRD_pairs_shared=trans_crd_cell2_CRD_pairs[(trans_crd_cell2_CRD_pairs$V4 %in% cell2shared) & (trans_crd_cell2_CRD_pairs$V8 %in% cell2shared), ]
  
  file.create(filename)
  # replace in trans_crd_cell1_CRD_pairs_shared with the equivalent in cell2 -B
  for (i in 1:length(trans_crd_cell1_CRD_pairs_shared[,1])) {
    crd_pair_cell1=c(trans_crd_cell1_CRD_pairs_shared[i,]$V4,trans_crd_cell1_CRD_pairs_shared[i,]$V8)
    cell2_equivalent_crd1 = shared_crds$V2[which(shared_crds$V1 == trans_crd_cell1_CRD_pairs_shared[i,]$V4)]
    cell2_equivalent_crd2 = shared_crds$V2[which(shared_crds$V1 == trans_crd_cell1_CRD_pairs_shared[i,]$V8)]
    crd_pair_cell1_in_cell2=c(cell2_equivalent_crd1 ,cell2_equivalent_crd2) # equivalent of cell1_pair in cell2_pair
    # pair_test=c("10_internal_7905",  "19_internal_6283") # idx  5160 : 10_internal_7905  19_internal_6283, from trans_crd_cell2_CRD_pairs_shared
    # find cell2_pair in trans_crd_cell2_CRD_pairs_shared
    indexes1=which(trans_crd_cell2_CRD_pairs_shared$V8 == crd_pair_cell1_in_cell2[1] &  trans_crd_cell2_CRD_pairs_shared$V4 == crd_pair_cell1_in_cell2[2])
    indexes2=which(trans_crd_cell2_CRD_pairs_shared$V4 == crd_pair_cell1_in_cell2[1] &  trans_crd_cell2_CRD_pairs_shared$V8 == crd_pair_cell1_in_cell2[2])
    for (index in c(indexes1,indexes2)){
      if (length(trans_crd_cell2_CRD_pairs_shared[index,1]) != 0){
        line=paste0(trans_crd_cell2_CRD_pairs_shared[index,]$V4,' ',trans_crd_cell2_CRD_pairs_shared[index,]$V8)
        write(line,file=filename,append=TRUE)
      }
    }
  }
  result=as.data.frame(data.table::fread(filename), head=FALSE, stringsAsFactors=FALSE)
  result 
}


#############################################################################################
#
# FUNCTION: PLOT
#
#############################################################################################

# without color limit (no max at 1)
plot_correlation_matrix_CRD_sharing_no_color_limit <- function(overlap_array_input, name, plot_directory){
  
  neut_vs_mono=overlap_array_input[1]
  neut_vs_tcell=overlap_array_input[2]
  mono_vs_neut=overlap_array_input[3]
  mono_vs_tcell=overlap_array_input[4]
  tcell_vs_mono=overlap_array_input[5]
  tcell_vs_neut=overlap_array_input[6]
  
  pdf(paste0(plot_directory,name,"_pairwise_comparisons.pdf"))
  M = matrix(c(1,neut_vs_mono,neut_vs_tcell,
               mono_vs_neut,1,mono_vs_tcell,
               tcell_vs_neut,tcell_vs_mono,1),
             ncol=3,byrow=T)
  colnames(M) = c("Neutrophils","Monocytes","T cells")
  rownames(M) = c("Neutrophils","Monocytes","T cells")
  corrplot(M, method = "number",is.corr=F,col = "black",number.cex=1.5)
  corrplot(M,is.corr=F,p.mat = M,sig.level=-1,insig = "p-value",number.cex=1.5)
  dev.off()
  
}


plot_correlation_matrix_CRD_sharing <- function(overlap_array_input, name, plot_directory, digits=4){
  
  neut_vs_mono=overlap_array_input[1]
  neut_vs_tcell=overlap_array_input[2]
  mono_vs_neut=overlap_array_input[3]
  mono_vs_tcell=overlap_array_input[4]
  tcell_vs_mono=overlap_array_input[5]
  tcell_vs_neut=overlap_array_input[6]
  
  pdf(paste0(plot_directory,name,"_pairwise_comparisons.pdf"))
  M = matrix(c(1,neut_vs_mono,neut_vs_tcell,
               mono_vs_neut,1,mono_vs_tcell,
               tcell_vs_neut,tcell_vs_mono,1),
             ncol=3,byrow=T)
  colnames(M) = c("Neutrophils","Monocytes","T cells")
  rownames(M) = c("Neutrophils","Monocytes","T cells")
  # corrplot(M,is.corr=F,cl.lim = c(0, 1),p.mat = M,sig.level=-1,insig = "p-value",number.cex=1.5)
  corrplot(M, method = "number",is.corr=F,col = "black",number.cex=1.5,number.digits=digits)
  dev.off()

}


#############################################################################################
#
# Folders and Files
#
#############################################################################################

path_shared_crds='/Users/dianaavalos/Programming/A_CRD_plots/CRD_sharing/'

path_transCRDs_signif = "/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant/"
# bis has start and end of CRDs
path_transCRDs_all = "/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:merged/"
path_out="/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:shared/"
plot_directory="/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:shared/plots/"

#############################################################################################
#
# MAIN 2
#
#############################################################################################

cell_pairs=c(c('neut','mono'),c('neut','tcell'),c('mono','neut'),c('mono','tcell'),c('tcell','neut'),c('tcell','mono'))
FDRthreshold=0.01

module='mean'
data_type='hist'
i=1
for(data_type in c('hist','methyl')){
  for(module in c('mean','loom')){
    for (i in c(1,3,5,7,9,11)){

      cell1=cell_pairs[i]
      cell2=cell_pairs[i+1]
      
      name=paste0(data_type,'_',module,'_',cell1,'_vs_',cell2)
      cat (name, ' ')
      
      trans_crd_cell1_signif=as.data.frame(data.table::fread(paste0(path_transCRDs_signif,data_type,'_', cell1,'_',module, '_trans.significant_',FDRthreshold,'.txt'), head=TRUE, stringsAsFactors=FALSE))
      trans_crd_cell2_signif=as.data.frame(data.table::fread(paste0(path_transCRDs_signif,data_type,'_', cell2,'_',module, '_trans.significant_',FDRthreshold,'.txt'), head=TRUE, stringsAsFactors=FALSE))
      trans_crd_cell1_all=as.data.frame(data.table::fread(paste0(path_transCRDs_all,data_type,'_', cell1,'_',module, '_trans_ALL.txt.gz'), head=TRUE, stringsAsFactors=FALSE))
      trans_crd_cell2_all=as.data.frame(data.table::fread(paste0(path_transCRDs_all,data_type,'_', cell2,'_',module, '_trans_ALL.txt.gz'), head=TRUE, stringsAsFactors=FALSE))
      shared_crds=as.data.frame(data.table::fread(paste0(path_shared_crds,name,'_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE))
      
      colnames(trans_crd_cell1_all)=colnames(trans_crd_cell2_all)=colnames(trans_crd_cell1_signif[,1:10])
      
      compute_shared_transCRD_v2(shared_crds,trans_crd_cell1_signif,trans_crd_cell2_all,name)

      
    }
  }
}

transCRD_c1_signif_shared=compute_shared_proportion_transCRD_assoc(trans_crd_cell1_signif,shared_crds$V1)
transCRD_c2_signif_shared=compute_shared_proportion_transCRD_assoc(trans_crd_cell2_signif,shared_crds$V2)
transCRD_c1_all_shared=compute_shared_proportion_transCRD_assoc(trans_crd_cell1_all,shared_crds$V1)
transCRD_c2_all_shared=compute_shared_proportion_transCRD_assoc(trans_crd_cell2_all,shared_crds$V2)

compute_shared_proportion_transCRD_assoc <- function(trans_crd_cell1_signif,cell1shared){
  # sharedCRD$V1 correspond to the CRDs of cell 1
  # cell1shared=shared_crds$V1
  trans_crd_cell1_CRD_pairs_shared=trans_crd_cell1_signif[(trans_crd_cell1_signif$idx1 %in% cell1shared) & (trans_crd_cell1_signif$idx2 %in% cell1shared), ]
  trans_crd_cell1_CRD_pairs_shared
}


#############################################################################################
#
# MAIN
#
#############################################################################################

cell_pairs=c(c('neut','mono'),c('neut','tcell'),c('mono','neut'),c('mono','tcell'),c('tcell','neut'),c('tcell','mono'))


for(FDRthreshold in c(0.01,0.05)){
  for(data_type in c('hist','methyl')){
    for(condition in c('mean','loom')){ 
      for (i in c(1,3,5,7,9,11)){
        cell1=cell_pairs[i]
        cell2=cell_pairs[i+1]
        cat(data_type, '  ',condition, '   ',cell1, '  ', cell2 )
        # load data
        trans_crd_cell1=as.data.frame(data.table::fread(paste0(path_transCRDs_signif,data_type,'_', cell1,'_',condition, '_trans.significant_',FDRthreshold,'.txt'), head=FALSE, stringsAsFactors=FALSE))
        trans_crd_cell2=as.data.frame(data.table::fread(paste0(path_transCRDs_signif,data_type,'_', cell2,'_',condition, '_trans.significant_',FDRthreshold,'.txt'), head=FALSE, stringsAsFactors=FALSE))
        shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,data_type,'_mean_',cell1,'_vs_',cell2,'_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE))
        # compute shared transCRDs
        result=compute_shared_transCRD(shared_crds,trans_crd_cell1,trans_crd_cell2,paste0(path_out,data_type,'_',condition,'_',cell1,'_vs_',cell2,'_',FDRthreshold,'_transCRDs_shared.txt'))
      }
    }
  }
}








#############################################################################################
#
# FURTHER ANALYSIS
#
#############################################################################################
cell_pairs[1:2]
for (cell1 in c('neut','mono','tcell')){
  for (cell2 in c('neut','mono','tcell')){
    if (cell1 != cell2){
    }}}

# Once that we have that we want to plot a corr matrix with the percentage of shared, and also just the numbers

for(FDRthreshold in c(0.01,0.05)){
  for(data_type in c('hist','methyl')){
    for(condition in c('mean','loom')){ 
      
      ### only the number of overlaps
      neut_vs_mono=length(readLines(paste0(path_out,data_type,'_',condition,'_','neut','_vs_','mono','_',FDRthreshold,'_transCRDs_shared.txt')))
      neut_vs_tcell=length(readLines(paste0(path_out,data_type,'_',condition,'_','neut','_vs_','tcell','_',FDRthreshold,'_transCRDs_shared.txt')))
      mono_vs_neut=length(readLines(paste0(path_out,data_type,'_',condition,'_','mono','_vs_','neut','_',FDRthreshold,'_transCRDs_shared.txt')))
      mono_vs_tcell=length(readLines(paste0(path_out,data_type,'_',condition,'_','mono','_vs_','tcell','_',FDRthreshold,'_transCRDs_shared.txt')))
      tcell_vs_mono=length(readLines(paste0(path_out,data_type,'_',condition,'_','tcell','_vs_','mono','_',FDRthreshold,'_transCRDs_shared.txt')))
      tcell_vs_neut=length(readLines(paste0(path_out,data_type,'_',condition,'_','tcell','_vs_','neut','_',FDRthreshold,'_transCRDs_shared.txt')))
      overlap_array_input=c(neut_vs_mono, neut_vs_tcell, mono_vs_neut, mono_vs_tcell , tcell_vs_mono, tcell_vs_neut)
      name=paste0(data_type,'_',condition,'_',FDRthreshold)
      cat (name,'   ', overlap_array_input)
      plot_correlation_matrix_CRD_sharing_no_color_limit(overlap_array_input, name, plot_directory)
      
      ### computing ratio
      # query_vs_ref, we want length(query) as denominator. see 23R
      # hist_neut_loo_trans.significant_0.05
      neut_vs_mono_r=neut_vs_mono/length(readLines(paste0(path_transCRDs_signif,data_type,'_neut_',condition,'_trans.significant_',FDRthreshold,'.txt')))
      neut_vs_tcell_r=neut_vs_tcell/length(readLines(paste0(path_transCRDs_signif,data_type,'_neut_',condition,'_trans.significant_',FDRthreshold,'.txt')))
      mono_vs_neut_r=mono_vs_neut/length(readLines(paste0(path_transCRDs_signif,data_type,'_mono_',condition,'_trans.significant_',FDRthreshold,'.txt')))
      mono_vs_tcell_r=mono_vs_tcell/length(readLines(paste0(path_transCRDs_signif,data_type,'_mono_',condition,'_trans.significant_',FDRthreshold,'.txt')))
      tcell_vs_mono_r=tcell_vs_mono/length(readLines(paste0(path_transCRDs_signif,data_type,'_tcell_',condition,'_trans.significant_',FDRthreshold,'.txt')))
      tcell_vs_neut_r=tcell_vs_neut/length(readLines(paste0(path_transCRDs_signif,data_type,'_tcell_',condition,'_trans.significant_',FDRthreshold,'.txt')))
      overlap_array_input_r=c(neut_vs_mono_r, neut_vs_tcell_r, mono_vs_neut_r, mono_vs_tcell_r , tcell_vs_mono_r, tcell_vs_neut_r)
      overlap_array_input_r_percent=overlap_array_input_r*100
      # in the others we did fraction, do we want percent here?
      name=paste0('ratio_',data_type,'_',condition,'_',FDRthreshold)
      cat (name,'   ', overlap_array_input_r)
      plot_correlation_matrix_CRD_sharing(overlap_array_input_r, name, plot_directory, digits=6)
    }
  }
}


