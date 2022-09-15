# compare the correlation based on independent groups  ------------------------------------
# create subfiles with sections to limit the ram

# cocor.indep.groups(r1.jk, r2.hm, n1, n2, alternative = "two.sided",
#                    test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
#                    data.name = NULL, var.labels = NULL, return.htest = FALSE)

# Clean environment ------------------------------------

rm(list = ls())
gc()


# Packages ---------------------------------------------

library(cocor)
library(dplyr)
library(data.table) # call fread
library(R.utils) # to read .gz

args = commandArgs(trailingOnly=TRUE)

# Paths ---------------------------------------------
#chr_num=args[1]

# path mac
path = '/Users/dianaavalos/Programming/A_CRD_plots/1_CRD'
path_peak = '/Users/dianaavalos/Programming/A_CRD_plots/8_PEAKS'
outpath= '/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS'


# path baobab to be updated
path = '/home/users/a/avalosma/scratch/1_CRD'
path_peak = '/home/users/a/avalosma/scratch/8_PEAKS'
outpath= '/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS'
path_script = "/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"

# Load files ---------------------------------------------

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94
sub <- 1
chr_num <- 5

for (chr_num in seq(1,22)){
  
  print(chr_num)
  
  CORR_neut = fread(file.path(path_peak, paste0(type, "_",  "neut", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE, select = c(1,2,3,7))
  CORR_mono = fread(file.path(path_peak, paste0(type, "_",  "mono", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE, select = c(1,2,3,7))
  CORR_tcell = fread(file.path(path_peak, paste0(type, "_",  "tcell", ".corr.chr", chr_num, ".txt")),
                     head = FALSE, stringsAsFactors = FALSE, select = c(1,2,3,7))
  
  
  # Main ---------------------------------------------
  
  # limit in cis the corr
  CORR_neut1 <- subset(CORR_neut, CORR_neut$V2 - CORR_neut$V1 <= 250)
  CORR_mono1 <- subset(CORR_mono, CORR_mono$V2 - CORR_mono$V1 <= 250)
  CORR_tcell1 <- subset(CORR_tcell, CORR_tcell$V2 - CORR_tcell$V1 <= 250)
  

  #  split the text in big chunks and process.
  # split the CORR files in 100 parts
  
  for (sub in seq(1,1000)){
    
    if (sub !=1000){
      index_1 <- round(dim(CORR_tcell1)[1]/100)*(sub-1)+1
      index_2 <- round(dim(CORR_tcell1)[1]/100)*sub
    } else{
      index_1 <- round(dim(CORR_tcell1)[1]/100)*(sub-1)+1
      index_2 <- dim(CORR_tcell1)[1]
    }
    
    CORR_neut_sub <- CORR_neut1[index_1:index_2,]
    CORR_mono_sub <- CORR_mono1[index_1:index_2,]
    CORR_tcell_sub <- CORR_tcell1[index_1:index_2,]
    
    CORR_neut_sub_file <- file.path(outpath, paste0("neut_chr",chr_num,"_section_",sub,".RData"))
    CORR_mono_sub_file <- file.path(outpath, paste0("mono_chr",chr_num,"_section_",sub,".RData"))
    CORR_tcell_sub_file <- file.path(outpath, paste0("tcell_chr",chr_num,"_section_",sub,".RData"))
    
    # save sub files
    save(CORR_neut_sub, file = CORR_neut_sub_file)
    save(CORR_mono_sub, file = CORR_mono_sub_file)
    save(CORR_tcell_sub, file = CORR_tcell_sub_file)

    }
}    






