# GOAL: here we want to plot the TFBS probability with the pvalue of the peak pair with -log10(pvalue) in X, and frequency (not count) in Y.
# we fornateed the peak dataframe and the TFBS dataframe in the functions here

# difference with v1: we have df4 <- fread("/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files/CORR_mono_neut_tcell_chr_22.txt")
# that contains cell 1 2 and 3 directly

# Inputs ---------------------------------------------
args = commandArgs(trailingOnly=TRUE)
chr=as.integer(args[1])
print(chr)
threshold=as.integer(args[2])

# Packages ---------------------------------------------

library(data.table) # for data.frame
library(tidyverse) # for <- %>%
library(mltools) # for bins

# Paths ---------------------------------------------

# mac
folder_CORR="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files"
outfolder_TF_pval_bins="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
peak_dir="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files/peaks_threshold"
folder_bins=outfolder_TF_pval_bins

# baobab
folder_CORR="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/files/merged"
peak_dir="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/peaks_threshold"

outfolder_TF_pval_bins="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
folder_plots="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/PLOTS"
folder_bins="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/bins"

# Load data ---------------------------------------------

name_file_intersect=paste0(peak_dir,"/intersect_chr_",chr,"_threshold_",threshold,".txt")
print(name_file_intersect)
bins_sec= c(10E-10,10E-5,10E-4,10E-3,10E-2,10E-1)


if (file.info(name_file_intersect)$size>1){
  
  file_intersect=fread(paste0(peak_dir,"/intersect_chr_",chr,"_threshold_",threshold,".txt"), head=TRUE)
  colnames(file_intersect) <- c("chr","startpeak","endpeak","peak","chrbis","startTF","endTF","origin","name")
  
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
    write.table(dout, file = paste0(folder_bins,"/table_",pair,"_",chr,".txt"), sep = "\t",row.names = FALSE, col.names = TRUE, quote=FALSE)
    
    d2 <- df %>% filter(pval_cell_1< 0.05)
    d3 <- df %>% filter(pval_cell_1< 0.05) %>% filter(pval_cell_2< 0.05)
    dfList <- list(d2,d3)
    
    # we are doing this over df3. maybe  try on df1 and df2
    index=1
    for (index in seq(1,length(dfList))){
      df_loop <- dfList[[index]]
      print(paste0("index ",index))  
      # count TF types per bin
      TFBS_unique_names <- data.frame(unique(file_intersect$name))
      colnames(TFBS_unique_names)="TF"
      
      #BINS LOOP
      for (j in seq(1,length(bins))){
        bin <- bins[j]
        bin_name <- paste0("[",bins_sec[j],",",bins_sec[j+1],"]")
        print(paste0(j," ",bin))
        dbin <- df_loop %>% filter(binned_pval_cell_1_2==bin) # here we are taking d3, this can be changed)
        peaks_in_bin <- unique( unlist(list(dbin$peak_i,dbin$peak_j)))
        TF_from_peaks <- file_intersect %>% filter(peak %in% peaks_in_bin) %>%  pull("name")
        if (length(TF_from_peaks)!=0){
          TF_count_bin <- table(TF_from_peaks) # how many of each TF in the bin
          TF_count_bin <- data.frame(TF_count_bin)
          colnames(TF_count_bin)=c("TF",toString(bin))
          number_of_TFs_detected <- sum(TF_count_bin[,2])
          TF_count_bin[,2] <- TF_count_bin[,2]/number_of_TFs_detected
          colnames(TF_count_bin) <- c("TF",bin_name)
          TFBS_unique_names <- merge(TFBS_unique_names,TF_count_bin, by.y ="TF",all = TRUE)
          TFBS_unique_names[is.na(TFBS_unique_names)] <- 0
        } else {
          TFBS_unique_names[bin_name] = 0.00
        }
        
        #head(TFBS_unique_names,10)
        write.table(TFBS_unique_names, file = paste0(outfolder_TF_pval_bins,"/TFs_d_",index,"_",pair,"_",chr,"_thre_",threshold,".txt"), sep = "\t",row.names = FALSE, col.names = TRUE, quote=FALSE)
        print(paste0("/TFs_d_",index,"_",pair,"_",chr,"_thre_",threshold,".txt"))
      }
    }
  }
} else {
  print ("empty")
}


