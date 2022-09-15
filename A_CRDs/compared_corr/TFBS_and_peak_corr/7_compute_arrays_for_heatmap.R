# Build heatmap for TFBS
# take the intersect file threshold 10 d2 per pair
# filter the TFs most represented
# build heatmap

# Params ---------------------------------------------
chr=22
pair="mono_neut"
peak_pair=1

threshold=0.01
threshold1=10


args = commandArgs(trailingOnly=TRUE)
chr=as.integer(args[1])
pair = args[2]

# Packages ---------------------------------------------

library(ggplot2)
library(data.table) # for data.frame
library(tidyverse) # for <- %>%
library(mltools) # for bins


# Load data ---------------------------------------------

# # mac
folder_CORR="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files"
peak_dir="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files/peaks_threshold" # intersect
TF_folder="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
outpath="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TF_arrays/"

# # baobab
folder_CORR="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/files/merged"
peak_dir="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/peaks_threshold"
TF_folder="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
outpath="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TF_arrays"
TFBSpath="/home/users/a/avalosma/scratch/Annotations/peaks_TFBS_overlap"

# 
# outfolder_TF_pval_bins="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
# folder_plots="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/PLOTS"
# folder_bins="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/bins"


# LOAD FILE INTERSECT ---------------------------------------------

name_file_intersect=paste0(peak_dir,"/intersect_chr_",chr,"_threshold_",threshold1,".txt")
print(name_file_intersect)
file_intersect=fread(name_file_intersect, head=FALSE)
colnames(file_intersect) <- c("chr","startpeak","endpeak","peak","chrbis","startTF","endTF","origin","name")

# LOAD FILE TFBS ---------------------------------------------

TFBS_peak = fread(file.path(TFBSpath,"TFBS_peaks_formatted.txt"), header=TRUE) # has all chromosomes already
length(unique(TFBS_peak$new_names)) # how many TFBS types ?
TFBS_ordered <- fread(paste0(TF_folder, "/TFs_most_represented_only_TF_ordered.txt"),header=FALSE)[1:50]
TFBS_ordered=TFBS_ordered$V1

# filter intersect with 50 most important

file_intersect <- file_intersect %>% filter(name %in% TFBS_ordered)

# MAIN ---------------------------------------------

pairs=c("mono_neut","mono_tcell","neut_tcell")
pairs=c("neut_mono","tcell_mono","tcell_neut")

if (pair=="mono_neut"){
  dout= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,7,8,11,12))
} else if (pair=="mono_tcell"){
  dout= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,9,10,13,14))
} else if (pair=="neut_tcell"){
  dout= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,7,8,9,10,15,16))
} else if (pair=="neut_mono"){
  dout= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,7,8,5,6,11,12))
} else if (pair=="tcell_mono"){
  dout= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,9,10,5,6,13,14))
} else {
  dout= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,9,10,7,8,15,16))
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

peak_pair=1
for (peak_pair in seq(1,dim(peaks_correlated_cell1)[1])){
  # find the TFBS associated with both peaks of the row and compute all the combinations
  TFBS_peak_i <- file_intersect %>% filter (peak %in% peaks_correlated_cell1[peak_pair,]$peak_i)
  TFBS_peak_j <- file_intersect %>% filter (peak %in% peaks_correlated_cell1[peak_pair,]$peak_j)
  
  # if both peaks overlap with TFBS meaning TFBS_peak_i and TFBS_peak_j not empty
  if (!dim(TFBS_peak_j)[1]==0 & !dim(TFBS_peak_i)[1]==0 ){
    # drop NA, which corresponds to TFBS that are not in the top 50, redundant because intersect file filtered?
    # replace names of TF_combinations by index numbers
    TFBS_peak_index_i <- discard(match(TFBS_peak_i$name,TFBS_ordered),is.na)
    TFBS_peak_index_j <- discard(match(TFBS_peak_j$name,TFBS_ordered),is.na)
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

print("finished")
