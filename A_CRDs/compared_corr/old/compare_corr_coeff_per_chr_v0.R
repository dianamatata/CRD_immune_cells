#!/usr/bin/env Rscript
# problem: files too big
# compare the correlation based on independent groups  ------------------------------------

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
chr_num=args[1]
print(chr_num)  

# path mac
path = '/Users/dianaavalos/Programming/A_CRD_plots/1_CRD'
path_peak = '/Users/dianaavalos/Programming/A_CRD_plots/8_PEAKS'

# path baobab to be updated
path = '/home/users/a/avalosma/scratch/1_CRD'
path_peak = '/home/users/a/avalosma/scratch/8_PEAKS'

# Functions  ---------------------------------

COR_c1_vs_c2 <- function(COR_c1,COR_c2,N1,N2) {
  COR_c1_vs_c2_df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(COR_c1_vs_c2_df) <- c("peak1", "peak2", "corr_c1","corr_c2", "fisher_statistic", "pvalue")
  for (i in seq(1,dim(COR_c1)[1])){
    a <- cocor.indep.groups(COR_c1$V7[i], COR_c2$V7[i], N1, N2, alternative = "two.sided",
                            test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                            data.name = NULL, var.labels = NULL, return.htest = FALSE)
    COR_c1_vs_c2_df <- rbind(COR_c1_vs_c2_df,c(COR_c1$V2[i],COR_c1$V3[i],COR_c1$V7[i],COR_c2$V7[i],a@fisher1925$statistic,a@fisher1925$p.value))
  }
  colnames(COR_c1_vs_c2_df) <- c("peak1", "peak2", "corr_c1","corr_c2", "fisher_statistic", "pvalue")
  COR_c1_vs_c2_df
}

# Load files ---------------------------------------------

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94

  
CORR_neut = fread(file.path(path_peak, paste0(type, "_",  "neut", ".corr.chr", chr_num, ".txt")),
                  head = FALSE, stringsAsFactors = FALSE)
CORR_mono = fread(file.path(path_peak, paste0(type, "_",  "mono", ".corr.chr", chr_num, ".txt")),
                  head = FALSE, stringsAsFactors = FALSE)
CORR_tcell = fread(file.path(path_peak, paste0(type, "_",  "tcell", ".corr.chr", chr_num, ".txt")),
                   head = FALSE, stringsAsFactors = FALSE)


# Main ---------------------------------------------

# limit in cis the corr
CORR_neut <- subset(CORR_neut, CORR_neut$V2 - CORR_neut$V1 <= 250)
CORR_mono <- subset(CORR_mono, CORR_mono$V2 - CORR_mono$V1 <= 250)
CORR_tcell <- subset(CORR_tcell, CORR_tcell$V2 - CORR_tcell$V1 <= 250)

# compute the compared correlation coeff for all combinaison    ---------------------------------------------
print("computing")


COR_mono_vs_neut <- COR_c1_vs_c2(CORR_mono, CORR_neut, Nmono, Nneut)
write.table(COR_mono_vs_neut, file =file.path(path, paste0("COR_mono_vs_neut_chr",chr_num,".txt")), sep = "\t",
            row.names = TRUE, col.names = NA)
remove(COR_mono_vs_neut)

COR_mono_vs_tcell <- COR_c1_vs_c2(CORR_mono, CORR_tcell, Nmono, Ntcell)
write.table(COR_mono_vs_tcell, file =file.path(path, paste0("COR_mono_vs_tcell_chr",chr_num,".txt")), sep = "\t",
            row.names = TRUE, col.names = NA)
remove(COR_mono_vs_tcell)

COR_neut_vs_tcell <- COR_c1_vs_c2(CORR_neut, CORR_tcell, Nneut, Ntcell)
write.table(COR_neut_vs_tcell, file =file.path(path, paste0("COR_neut_vs_tcell_chr",chr_num,".txt")), sep = "\t",
            row.names = TRUE, col.names = NA)
remove(COR_neut_vs_tcell)
  
print("finished")
