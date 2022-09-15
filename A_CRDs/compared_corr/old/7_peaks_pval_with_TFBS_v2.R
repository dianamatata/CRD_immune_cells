# GOAL: here we want to plot the TFBS probability with the pvalue of the peak pair with -log10(pvalue) in X, and frequency (not count) in Y.
# we fornateed the peak dataframe and the TFBS dataframe in the functions here

# difference with v1: we have df4 <- fread("/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files/CORR_mono_neut_tcell_chr_22.txt")
# that contains cell 1 2 and 3 directly

# Inputs ---------------------------------------------
args = commandArgs(trailingOnly=TRUE)
chr=as.integer(args[1])
print(chr)
# Packages ---------------------------------------------

library(data.table) # for data.frame
library(tidyverse) # for <- %>% 
library(mltools) # for bins

# Paths ---------------------------------------------

# mac
folder_TFBS="/Users/dianaavalos/Desktop/TFBS"
folder_CORR="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files"
folder_TF_pval_bins="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TFBS"

# baobab
folder_TFBS="/home/users/a/avalosma/scratch/Annotations/peaks_TFBS_overlap"
folder_CORR="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/files/merged"
folder_TF_pval_bins="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
folder_plots="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/PLOTS"
folder_bins="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/bins"

# Load data ---------------------------------------------
TFBS_peak = fread(file.path(folder_TFBS,"TFBS_peaks_formatted.txt"), header=TRUE) # has all chromosomes already
length(unique(TFBS_peak$new_names)) # how many TFBS types ?

pairs=c("mono_neut","mono_tcell","neut_tcell")
for (pair in pairs){
    print(paste0("processing ",pair,"_",chr))
    
    if (pair=="mono_neut"){
      df= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,7,8,11,12))
    } else if (pair=="mono_tcell"){
      df= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,5,6,9,10,13,14))
    } else {
      df= fread(paste0(folder_CORR,"/CORR_mono_neut_tcell_chr_",chr,"_ALL_sorted.txt"), header=TRUE, select=c(1,2,3,4,7,8,9,10,15,16))
    }
    colnames(df) <- c("peak_i","peak_j","location_i","location_j","rho_cell_1","rho_cell_2","pval_cell_1","pval_cell_2","rho_cell_1_2","pval_cell_1_2")
    
    # create bins and count per bin
    df$binned_pval_cell_1_2 <- bin_data(df$pval_cell_1_2, bins= c(10E-10,10E-5,10E-4,10E-3,10E-2,10E-1), boundaryType = "[lorc") # adding bins in pvalue
    bins <- sort(unique(df$binned_pval_cell_1_2 ))
    dout <- table(df$binned_pval_cell_1_2)
    write.table(df, file = paste0(folder_bins,"/",pair,"_",chr,".txt"), sep = "\t",row.names = FALSE, col.names = TRUE, quote=FALSE)
    d2 <- df %>% filter(pval_cell_1< 0.05)
    d3 <- df %>% filter(pval_cell_1< 0.05) %>% filter(pval_cell_2< 0.05)
    
    # we are doing this over df3. maybe later try on df1 and df2
    
    # count TF types per bin
    TFBS_unique_names <- data.frame(unique(TFBS_peak$new_names))
    colnames(TFBS_unique_names)="TF"
    for (j in seq(1,length(bins))){
      print(j)
      bin=bins[j]
      print(bin)
      dbin <- d3 %>% filter(binned_pval_cell_1_2==bin) # here we are taking d3, this can be changed
      peaks_in_bin <- unique( unlist(list(dbin$peak_i,dbin$peak_j))) 
      TFBS_peak_nbr <- TFBS_peak %>% filter (peak_nbr %in% peaks_in_bin)
      TF_count_bin <- table(TFBS_peak_nbr$new_names)
      TF_count_bin <- data.frame(TF_count_bin)
      colnames(TF_count_bin)=c("TF",toString(bin))
      number_of_TFs_detected <- sum(TF_count_bin[,2])
      # TF_count_bin$freq <- TF_count_bin[,2]/number_of_TFs_detected
      TF_count_bin[,2] <- TF_count_bin[,2]/number_of_TFs_detected
      head(TF_count_bin)
      TFBS_unique_names <- merge(TFBS_unique_names,TF_count_bin, by.y ="TF",all = TRUE)
      TFBS_unique_names[is.na(TFBS_unique_names)] <- 0
    }
    #head(TFBS_unique_names,10)
    write.table(TFBS_unique_names, file = paste0(folder_TF_pval_bins,"/TFs_",pair,"_",chr,".txt"), sep = "\t",row.names = FALSE, col.names = TRUE, quote=FALSE)
    print(paste0("TFs_",pair,"_",chr," is written"))
}

