# args = commandArgs(trailingOnly=TRUE)
# chr_num=as.integer(args[1])
# threshold=as.integer(args[2])

# Packages ---------------------------------------------

library("dplyr")
library("data.table")
library("tidyverse")

# Paths ---------------------------------------------

datadir="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/files"
datadir="/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/files"
outdir="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/peaks_threshold"
# Main ---------------------------------------------

pairs=c("mono_neut","mono_tcell","neut_tcell")

chr=22
for (chr in seq(1,22)){
  file_peaks=fread(file.path(datadir,paste0("mono_neut_tcell_CIS.chr",chr,".txt")), head=FALSE)
  file_peaks$chr <- chr
  file_peaks$chr <- paste0('chr',file_peaks$chr)
  df <- unique(file_peaks[,c(11,3,1)])
  df <- df[order(V1),]
  # filter unique V3
  
  for (threshold in c(10,50,100,150)){
    df$start <- df$V3-threshold
    df$stop <- df$V3+threshold
    file_peaks_out <- df[,c(1,4,5,3)]
    colnames(file_peaks_out) <- c("chr" ,"start", "stop","peak")
    write.table(file_peaks_out, file = paste0(outdir,"/peaks_chr_",chr,"_threshold_",threshold,".txt"), sep = "\t",row.names = FALSE, quote = FALSE)
  }

}
