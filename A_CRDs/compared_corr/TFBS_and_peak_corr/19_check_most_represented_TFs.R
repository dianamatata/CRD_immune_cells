# open TFBS

library(ggplot2)
library(data.table) # for data.frame
library(tidyverse) # for <- %>%

# Paths ---------------------------------------------
data_annot="/Users/dianaavalos/Desktop/TFBS"

data_annot="/home/users/a/avalosma/scratch/Annotations"


# Main ---------------------------------------------
TFmotifmap=fread(file.path(data_annot,"HUMAN_hg19_BBLS_1_00_FDR_0_10_v1.bed.gz"), head=FALSE)
TFremap=fread(file.path(data_annot,"subset_remap2022_nr_macs2_hg19_v1_0.bed.gz"), head=FALSE)
TFremap=fread(file.path(data_annot,"remap2022_nr_macs2_hg19_v1_0.bed.gz"), head=FALSE)


TFmotifmap <- TFmotifmap[,1:4]
TFmotifmap$origin <- "motifmap"
TFmotifmap$new_names <- sapply(strsplit(TFmotifmap$V4, "="),function(x) x[[2]])

TFremap <- TFremap[,1:4]
TFremap$origin <- "remap"
TFremap$new_names <- sapply(strsplit(TFremap$V4, ":"),function(x) x[[1]])

TFBS <- rbind(TFremap,TFmotifmap)

# grep TF present in
TF_folder="/Users/dianaavalos/Programming/A_CRD_plots/CRD_project_outputs/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
TF_folder="/srv/beegfs/scratch/users/a/avalosma/CRD_project_outputs/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
TFBS_ordered <- fread(paste0(TF_folder, "/TFs_most_represented_only_TF_ordered.txt"),header=FALSE)[1:100]
TFBS_ordered=TFBS_ordered$V1

a <- TFBS %>% filter(new_names %in% TFBS_ordered)
write.table(a, file = paste0(data_annot,"/TFs_highly_present_100.txt"), sep = "\t",row.names = FALSE, quote = FALSE)


