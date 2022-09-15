# carefull chr from 2 to 22 not 1

# merge 
library(data.table) # for data.frame
library(tidyverse) # for <- %>% 
folder="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TF_arrays"
pairs=c("mono_neut","mono_tcell","neut_tcell")

for (pair in pairs){
  
  for (chr in seq(1,22)){
    print(paste0(pair, " ", chr))
    df_NOT_specific_NO_overlap <- read.csv(paste0(folder,"/",pair,"_",chr,"_array_NOT_specific_NO_overlap.txt"), row.names = 1, header= TRUE, sep="\t")
    df_specific_NO_overlap <- read.csv(paste0(folder,"/",pair,"_",chr,"_array_specific_NO_overlap.txt"), row.names = 1, header= TRUE, sep="\t")
    df_NOT_specific <- read.csv(paste0(folder,"/",pair,"_",chr,"_array_NOT_specific.txt"), row.names = 1, header= TRUE, sep="\t")
    df_specific <- read.csv(paste0(folder,"/",pair,"_",chr,"_array_specific.txt"), row.names = 1, header= TRUE, sep="\t")
    
    if (chr==1){
      ALL_spec <- df_specific
      ALL_NOT_spec <- df_NOT_specific
      ALL_spec_NO_O <- df_specific_NO_overlap
      ALL_NOT_spec_NO_O <- df_NOT_specific_NO_overlap
    } else {
      ALL_spec <- ALL_spec + df_specific
      ALL_NOT_spec <- ALL_NOT_spec + df_NOT_specific
      ALL_spec_NO_O <- ALL_spec_NO_O + df_specific_NO_overlap
      ALL_NOT_spec_NO_O <- ALL_NOT_spec_NO_O + df_NOT_specific_NO_overlap
    }
  }
  
  write.table(ALL_spec,
    file = paste0(folder, "/pairs/","ALL_spec_", pair, "_ALL.txt"),
    sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  write.table(ALL_spec_NO_O,
              file = paste0(folder, "/pairs/","ALL_spec_NO_O_", pair, "_ALL.txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  write.table(ALL_NOT_spec,
              file = paste0(folder, "/pairs/","ALL_NOT_spec_", pair, "_ALL.txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  write.table(ALL_NOT_spec_NO_O,
              file = paste0(folder, "/pairs/","ALL_NOT_spec_NO_O_", pair, "_ALL.txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}



