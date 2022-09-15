# Goal: add start end to TRH significant files
# so far  this one

# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(data.table)
library(tidyverse)


#############################################################################################
#
# DIRECTORIES AND FILES
#
#############################################################################################

directory='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant/'
crd_info_dir='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/CRD_id_start_end/'
out_dir='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant_bis/'


#############################################################################################
#
# MAIN
#
#############################################################################################
# methyl_neut_CRD_id.txt
# problem so far: same value for all starts and ends unique(TRH_signif$start1) only one value
files <- list.files(path=directory, pattern="0.0*.txt", full.names=TRUE, recursive=FALSE)


for(cell in c('neut','mono','tcell')){
  for(data_type in c('hist','methyl')){
    for(module in c('mean','loom')){
      for(FDR in c('0.01','0.05')){
        
        #file names
        name_CRD=paste0(data_type,'_',cell,'_CRD_id.txt')
        name_trans_file=paste0(data_type,'_',cell,'_',module,'_trans.significant_',FDR,'.txt')
        cat(name_trans_file,'  ',name_CRD,'  ')
        
        # load files
        CRDinfo = as.data.frame(data.table::fread(paste0(crd_info_dir,name_CRD), head=TRUE, stringsAsFactors=FALSE))
        TRH_signif_temp = as.data.frame(data.table::fread(paste0(directory,name_trans_file), head=TRUE, stringsAsFactors=FALSE))
        
        # reorder columns by name in  TRH_signif_temp
        TRH_signif_temp$start1=TRH_signif_temp$start2=TRH_signif_temp$end1=TRH_signif_temp$end2=NA
        col_order <- c("idx1","chr1","start1","end1","id1","idx2","chr2","start2","end2","id2","corr","pval","qval","midplace","midplace2")
        TRH_signif <- TRH_signif_temp[col_order]
        
        # update loop with start and stop
        
        for (i in 1:length(TRH_signif[,1])){
          idx_1=which(CRDinfo$id == TRH_signif[i,]$id1)
          TRH_signif$start1[i]=CRDinfo$start[idx_1]
          TRH_signif$end1[i]=CRDinfo$end[idx_1]
          idx_2=which(CRDinfo$id == TRH_signif[i,]$id2)
          TRH_signif$start2[i]=CRDinfo$start[idx_2]
          TRH_signif$end2[i]=CRDinfo$end[idx_2]
        }
        
        write.table(TRH_signif, paste0(out_dir,name_trans_file), append = FALSE,
                    sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
        
      }
    }
  }
}

# debug
name_trans_file='hist_neut_mean_trans.significant_0.01.txt' 
name_CRD='hist_neut_CRD_id.txt'
unique(TRH_signif$id1)
unique(TRH_signif$start1)

