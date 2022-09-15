# Packages ---------------------------------------------

library(ggplot2)
library(data.table) # for data.frame
library(tidyverse) # for <- %>%

# Paths ---------------------------------------------
data_annot="/Users/dianaavalos/Desktop/TFBS"

data_annot="/home/users/a/avalosma/scratch/Annotations"


# Main ---------------------------------------------
TFmotifmap=fread(file.path(data_annot,"HUMAN_hg19_BBLS_1_00_FDR_0_10_v1.bed"), head=FALSE)
TFremap=fread(file.path(data_annot,"remap2022_nr_macs2_hg19_v1_0.bed"), head=FALSE)


TFmotifmap <- TFmotifmap[,1:4]
TFmotifmap$origin <- "motifmap"
TFmotifmap$new_names <- sapply(strsplit(TFmotifmap$V4, "="),function(x) x[[2]])
TFmotifmap <- TFmotifmap[,-4]

TFremap <- TFremap[,1:4]
TFremap$origin <- "remap"
TFremap$new_names <- sapply(strsplit(TFremap$V4, ":"),function(x) x[[1]])
TFremap <- TFremap[,-4]

TFBS <- rbind(TFremap,TFmotifmap)
TFBS$V1 <- paste0('chr',TFBS$V1)

write.table(TFBS, file = paste0(data_annot,"/TFBS.txt"), sep = "\t",row.names = FALSE, quote = FALSE)


