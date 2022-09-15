# Build heatmap for TFBS
# The final idea is that the heatmap presents the odd ratio of each association of TFBS


# Packages ---------------------------------------------

library(ggplot2)
library(data.table) # for data.frame
library(tidyverse) # for <- %>% 
library(mltools) # for bins
# library('plot.matrix')


# Functions ---------------------------------------------


formatting_dataframe_2celltypes_correlation <- function(df){
  num_cols=12
  num_rows <- dim(df)[1]/num_cols
  dout=data.frame(matrix(ncol=num_cols,nrow=0))
  for(i in seq(0,num_rows-1)){    dout[dim(dout)[1]+1,] <- df$V1[(num_cols*i+1):(num_cols*(i+1))] }
  colnames(dout)=c("peak_i","peak_j","location_i","location_j","rho_cell_1","rho_cell_2","pval_cell_1","pval_cell_2","N1","N2","rho_cell_1_2","pval_cell_1_2")
  dout
}

# external arguments ---------------------------------------------
chr=22
pair="mono_neut"

args = commandArgs(trailingOnly=TRUE)
chr = args[1]
pair = args[2]


# Load data ---------------------------------------------


# mac
TFBSpath="/Users/dianaavalos/Desktop/TFBS"
folder_CORR="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files"
outpath="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TF_arrays"


# baobab
TFBSpath="/home/users/a/avalosma/scratch/Annotations/peaks_TFBS_overlap"
folder_CORR="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/files/merged"
outpath="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TF_arrays"


TFBS_peak = fread(file.path(TFBSpath,"TFBS_peaks_formatted.txt"), header=TRUE) # has all chromosomes already
length(unique(TFBS_peak$new_names)) # how many TFBS types ?


threshold=0.01
TFBS_ordered <- sort(unique(TFBS_peak$new_names))
peak_pair=1

#  TODO: TAKE ONLY TAKE 50 MOST REP
a <- data.frame(sort(table(TFBS_peak$new_names),decreasing=TRUE))
TFBS_ordered <- data.frame(sort(a[1:50,]$Var1))
colnames(TFBS_ordered) <- "TF"
TFBS_ordered <- TFBS_ordered$TF

# TFBS_ordered <- data.frame(sort(unique(TFBS_peak$new_names)))
# subset to train
# dout <- dout[1:610,]

pairs=c("mono_neut","mono_tcell","neut_tcell")
#for (pair in pairs){
 # for (chr in seq(1,22)){
    if (pair=="mono_neut"){
      df= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,7,8,11,12))
    } else if (pair=="mono_tcell"){
      df= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,9,10,13,14))
    } else {
      df= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,7,8,9,10,15,16))
    }
    colnames(df) <- c("peak_i","peak_j","location_i","location_j","rho_cell_1","rho_cell_2","pval_cell_1","pval_cell_2","rho_cell_1_2","pval_cell_1_2")
      
      #  Main  ---------------------------------------------
      
      # create sub arrays of peaks
      peaks_correlated_cell1 <- df %>% filter(pval_cell_1 <= threshold)
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
 # }
#}

