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
# chr_num=args[1]
chr_num=5
#
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

# Function ---------------------------------------------

COR_c1_vs_c2 <- function(COR_c1,COR_c2,N1,N2,filename) {
  header_txt<- c("peak1", "peak2", "corr_c1","corr_c2", "fisher_statistic", "pvalue","\n")
  cat(paste(unlist(header_txt),collapse="\t"),file=filename,sep="\n")
  
  for (i in seq(1,dim(COR_c1)[1])){
    a <- cocor.indep.groups(COR_c1$V7[i], COR_c2$V7[i], N1, N2, alternative = "two.sided",
                            test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                            data.name = NULL, var.labels = NULL, return.htest = FALSE)
    cat(paste(unlist(c(COR_c1$V2[i],COR_c1$V3[i],COR_c1$V7[i],COR_c2$V7[i],a@fisher1925$statistic,a@fisher1925$p.value)),collapse="\t"),file=filename,sep="\n",append=TRUE)
  }
}

# Load files ---------------------------------------------

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94
sub <- 1


CORR_neut = fread(file.path(path_peak, paste0(type, "_",  "neut", ".corr.chr", chr_num, ".txt")),
                  head = FALSE, stringsAsFactors = FALSE, select = c(1,2,3,7))
CORR_mono = fread(file.path(path_peak, paste0(type, "_",  "mono", ".corr.chr", chr_num, ".txt")),
                  head = FALSE, stringsAsFactors = FALSE, select = c(1,2,3,7))
CORR_tcell = fread(file.path(path_peak, paste0(type, "_",  "tcell", ".corr.chr", chr_num, ".txt")),
                   head = FALSE, stringsAsFactors = FALSE, select = c(1,2,3,7))


# Main ---------------------------------------------

# limit in cis the corr
CORR_neut <- subset(CORR_neut, CORR_neut$V2 - CORR_neut$V1 <= 250)
CORR_mono <- subset(CORR_mono, CORR_mono$V2 - CORR_mono$V1 <= 250)
CORR_tcell <- subset(CORR_tcell, CORR_tcell$V2 - CORR_tcell$V1 <= 250)

# 2 possibilities. write in file progressively with cat("World",file="outfile.txt",append=TRUE)

# test
CORR_neut <-  head(CORR_neut)
CORR_mono <-  head(CORR_mono)

# Possibility 1  ---------------------------------------------
filename=file.path(outpath, paste0("COR_mono_vs_neut_chr",chr_num,".txt"))
COR_c1_vs_c2(CORR_mono,CORR_neut,Nmono,Nneut,filename) 
  





  
