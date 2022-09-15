# after compare_corr_coeff in which we compute
# Clean environment ------------------------------------


# Packages ---------------------------------------------

library(cocor)
library(dplyr)
library(data.table) # call fread
library(R.utils) # to read .gz


# Paths ---------------------------------------------

# path mac
path = '/Users/dianaavalos/Programming/A_CRD_plots/1_CRD'
path_peak = '/Users/dianaavalos/Programming/A_CRD_plots/8_PEAKS'
path_corr = '/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/'
outpath = "/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/"

# path baobab to be updated
path = '/home/users/a/avalosma/scratch/1_CRD'
path_peak = '/home/users/a/avalosma/scratch/8_PEAKS'
path_corr = '/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/merged'
outpath = "/home/users/a/avalosma/scratch/8_PEAKS/PEAK_CIS"



# Main ---------------------------------------------

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94
chr_num <- 5



for (chr_num in 1:22){
  
  # Load files ---------------------------------------------
  CORR_neut = fread(file.path(path_peak, paste0(type, "_",  "neut", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  CORR_mono = fread(file.path(path_peak, paste0(type, "_",  "mono", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  CORR_tcell = fread(file.path(path_peak, paste0(type, "_",  "tcell", ".corr.chr", chr_num, ".txt")),
                     head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  
    # limit in cis the corr  ---------------------------------------------
  
  CORR_neut <- subset(CORR_neut, CORR_neut$V2 - CORR_neut$V1 <= 250)
  CORR_mono <- subset(CORR_mono, CORR_mono$V2 - CORR_mono$V1 <= 250)
  CORR_tcell <- subset(CORR_tcell, CORR_tcell$V2 - CORR_tcell$V1 <= 250)
  
  write.table(CORR_neut, file =file.path(outpath, paste0("CORR_neut",chr_num,".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  write.table(CORR_mono, file =file.path(outpath, paste0("CORR_mono",chr_num,".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  write.table(CORR_tcell, file =file.path(outpath, paste0("CORR_tcell",chr_num,".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
}

