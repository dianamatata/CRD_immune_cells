# take 50 most represented TFs
# subset these in each pair
# plot matplot with name of rows

library(data.table) # for data.frame
library(tidyverse) # for <- %>%
library(dplyr)

TF_folder="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
pairs=c("mono_neut","mono_tcell","neut_tcell")
TFs_to_keep <- fread(paste0(TF_folder, "/TFs_most_represented_only_TF_ordered.txt"),header=FALSE)[1:50]
for (pair in pairs){
  df_pair <- fread(file=paste0(TF_folder,"/TFs_",pair,"_ALL.txt"), header=TRUE)
  subset_df_pair <- df_pair %>% filter(TF %in% TFs_to_keep$V1)
  
  write.table(
    subset_df_pair,
    file = paste0(TF_folder, "/TFs_most_represented_",pair,"_v2.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

  
