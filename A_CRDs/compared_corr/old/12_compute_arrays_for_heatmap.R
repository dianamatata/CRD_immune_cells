# Build heatmap for TFBS
# The final idea is that the heatmap presents the odd ratio of each association of TFBS

# external arguments ---------------------------------------------
chr=22
pair="mono_neut"

args = commandArgs(trailingOnly=TRUE)
chr = args[1]
pair = args[2]

# Packages ---------------------------------------------

library(ggplot2)
library(data.table) # for data.frame
library(tidyverse) # for <- %>% 
library(mltools) # for bins
# library('plot.matrix')


# Functions ---------------------------------------------

formatting_TFBS_data <- function(TFBS_peak){
  g <- TFBS_peak %>% filter(origin == "motifmap") 
  g$new_names <- sapply(strsplit(g$TF_name, "="),function(x) x[[2]])
  h <- TFBS_peak %>% filter(origin == "remap") 
  h$new_names <- sapply(strsplit(h$TF_name, ":"),function(x) x[[1]])
  TFBS_peak <- rbind(h,g)
  TFBS_peak
}


# Load data ---------------------------------------------

# baobab
TF_folder="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
datadir="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/files/merged"
outpath="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TF_arrays"
TFBSpath="/home/users/a/avalosma/scratch/Annotations/peaks_TFBS_overlap"

TFBS_peak = fread(file.path(TFBSpath,"TFBS_peaks_formatted.txt"), header=TRUE) # has all chromosomes already
length(unique(TFBS_peak$new_names)) # how many TFBS types ?


TFBS_ordered <- fread(paste0(TF_folder, "/TFs_most_represented_only_TF_ordered.txt"),header=FALSE)[1:50]
TFBS_ordered=TFBS_ordered$V1
threshold=0.01
peak_pair=1


# TFBS_ordered <- data.frame(sort(unique(TFBS_peak$new_names)))
# subset to train
# dout <- dout[1:610,]

pairs=c("mono_neut","mono_tcell","neut_tcell")
#for (pair in pairs){
 #for (chr in seq(1,22)){
    if (pair=="mono_neut"){
      dout= fread(paste0(datadir,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,7,8,11,12))
    } else if (pair=="mono_tcell"){
      dout= fread(paste0(datadir,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,9,10,13,14))
    } else {
      dout= fread(paste0(datadir,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,7,8,9,10,15,16))
    }
    colnames(dout) <- c("peak_i","peak_j","location_i","location_j","rho_cell_1","rho_cell_2","pval_cell_1","pval_cell_2","rho_cell_1_2","pval_cell_1_2")
      
    #  Main  ---------------------------------------------
    
    # create sub arrays of peaks
    peaks_correlated_cell1 <- dout %>% filter(pval_cell_1 <= threshold)
    cell_specific_peaks_corr_in_c1 <- peaks_correlated_cell1%>% filter(pval_cell_1_2 <= threshold)
    not_cell_specific_peaks_corr_in_c1 <- peaks_correlated_cell1%>% filter(pval_cell_1_2 > threshold)
    
    print(paste0("size_peaks_cell_1_specific: ",dim(peaks_correlated_cell1)[1]))
    print(paste0("size_peak_pairs_specific_for_2: ",dim(cell_specific_peaks_corr_in_c1)[1]))
    print(paste0("size_peak_pairs_NOT_specific_for_2: ",dim(not_cell_specific_peaks_corr_in_c1)[1]))
    
    # create arrays to fill

    size_heatmap <- length(TFBS_ordered)
    array_specific = array_NOT_specific=matrix(0L,nrow=size_heatmap,ncol = size_heatmap)
    rownames(array_specific) <- TFBS_ordered
    colnames(array_specific) <- TFBS_ordered
    rownames(array_NOT_specific) <- TFBS_ordered
    colnames(array_NOT_specific) <- TFBS_ordered
    
    for (peak_pair in seq(1,dim(peaks_correlated_cell1)[1])){
      # find the TFBS associated with both peaks of the row and compute all the combinations
      # print(peaks_correlated_cell1[peak_pair,])
      TFBS_peak_i <- TFBS_peak %>% filter (peak_nbr %in% peaks_correlated_cell1[peak_pair,]$peak_i)
      TFBS_peak_j <- TFBS_peak %>% filter (peak_nbr %in% peaks_correlated_cell1[peak_pair,]$peak_j)
      
      # if both peaks overlap with TFBS meaning TFBS_peak_i and TFBS_peak_j not empty
      if (!dim(TFBS_peak_j)[1]==0 & !dim(TFBS_peak_i)[1]==0 ){
        # drop NA, which corresponds to TFBS that are not in the top 50
        # replace names of TF_combinations by index numbers
        TFBS_peak_index_i <- discard(match(TFBS_peak_i$new_names,TFBS_ordered),is.na)
        TFBS_peak_index_j <- discard(match(TFBS_peak_j$new_names,TFBS_ordered),is.na)
        TF_combinations <- unique(rbind(expand.grid(TFBS_peak_index_i,TFBS_peak_index_j),expand.grid(TFBS_peak_index_j,TFBS_peak_index_i))) # is symmetric
        
        # create a matrix for the TFBS
        current_matrix <- matrix(0L,nrow=size_heatmap,ncol = size_heatmap)
        for (row_i in seq(1,length(TF_combinations$Var1))){
          i <- TF_combinations[row_i,]$Var1
          j <- TF_combinations[row_i,]$Var2
          current_matrix[i,j] <- current_matrix[i,j] +1
        }
        # add this matrix to cell specific or not array
        if (peaks_correlated_cell1[peak_pair, ]$pval_cell_1_2 < threshold){
          array_specific <- array_specific + current_matrix       # significant
        } else{
          # NOT significant
          array_NOT_specific <- array_NOT_specific + current_matrix
        }
      }
    }
    
    # # how many peak pairs total are cell specific and not cell specific
    # dim(cell_specific_peaks_corr_in_c1)[1]
    # dim(not_cell_specific_peaks_corr_in_c1)[1]
    
    # then fill the not overlapping arrays
    array_specific_NO_overlap <- matrix(dim(cell_specific_peaks_corr_in_c1)[1],nrow=size_heatmap,ncol = size_heatmap) - array_specific
    array_NOT_specific_NO_overlap =matrix(dim(not_cell_specific_peaks_corr_in_c1)[1],nrow=size_heatmap,ncol = size_heatmap)  - array_NOT_specific
    
    rownames(array_specific_NO_overlap) = colnames(array_specific_NO_overlap) = rownames(array_NOT_specific_NO_overlap) = colnames(array_NOT_specific_NO_overlap) = TFBS_ordered

    # write files
    write.table(array_specific, file = paste0(outpath,"/",pair,"_",chr,"_array_specific.txt"), sep = "\t",row.names = TRUE, col.names = TRUE, quote=FALSE)
    write.table(array_NOT_specific, file = paste0(outpath,"/",pair,"_",chr,"_array_NOT_specific.txt"), sep = "\t",row.names = TRUE, col.names = TRUE, quote=FALSE)
    write.table(array_specific_NO_overlap, file = paste0(outpath,"/",pair,"_",chr,"_array_specific_NO_overlap.txt"), sep = "\t",row.names = TRUE, col.names = TRUE, quote=FALSE)
    write.table(array_NOT_specific_NO_overlap, file = paste0(outpath,"/",pair,"_",chr,"_array_NOT_specific_NO_overlap.txt"), sep = "\t",row.names = TRUE, col.names = TRUE, quote=FALSE)
  
 #}
#}


