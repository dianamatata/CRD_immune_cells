# compare the correlation based on independent groups  ------------------------------------

# cocor.indep.groups(r1.jk, r2.hm, n1, n2, alternative = "two.sided",
#                    test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
#                    data.name = NULL, var.labels = NULL, return.htest = FALSE)

# Packages ---------------------------------------------

library(cocor)
library(dplyr)
library(data.table) # call fread
library(R.utils) # to read .gz


# Paths ---------------------------------------------

# path mac
path = '/Users/dianaavalos/Programming/A_CRD_plots/1_CRD/'
outpath = "/Users/dianaavalos/Programming/A_CRD_plots/CRDs_with_diff_corr"


# path baobab to be updated
path = '/home/users/a/avalosma/scratch/1_CRD'
path_peak="/home/users/a/avalosma/scratch/8_PEAKS"
outpath='/home/users/a/avalosma/scratch/13_CRD_corr'

# Functions  ---------------------------------


MODS_only_CRDs <- function(MODS) {
  MODS1 = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ] # otherwise too many triangles
  MODS1
}


compute_corr_coeff_CRDs_and_significance <- function(COR_c1,COR_c2, MOD_cell_1) {
  
  summary <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(summary) <- c("chr_num", "CRD_name","CRD_start","CRD_end", "corr_c1","corr_c2", "fisher_statistic", "pvalue") 
  for (crd_i in seq(1, dim(MOD_cell_1)[1]) ) {
    peak_i=MOD_cell_1$LIDX[crd_i]
    peak_j=MOD_cell_1$RIDX[crd_i]
    
    # take all the correlation points between these 2 peaks for cell 1
    COR_cell_1 <- COR_c1 %>% filter(V1 >= peak_i & V2 <= peak_j)
    corr_coeff_mean_CRD_1 <- mean(COR_cell_1$V7)
    
    # do the same for cell type 2 but with CRDs of cell type 1
    COR_cell_2 <- COR_c2 %>% filter(V1 >= peak_i & V2 <= peak_j)
    corr_coeff_mean_CRD_2 <- mean(COR_cell_2$V7)
    
    N1=N2=dim(COR_cell_1)[1] # number of peaks in the CRD
    
    # compute correlation statistic between these 2 CRDs
    # what do we take as N? sample size or CRD size which is the same.
    a <- cocor.indep.groups(corr_coeff_mean_CRD_1, corr_coeff_mean_CRD_2, N1, N2, alternative = "two.sided",
                            test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                            data.name = NULL, var.labels = NULL, return.htest = FALSE)
    CRD_name=MOD_cell_1$UUID[crd_i]
    # maybe add peak start peak end to make the drawing if CRDs easier
    CRD_start <-  MOD_cell_1$LIDX[crd_i]
    CRD_end <-  MOD_cell_1$RIDX[crd_i]
    summary[dim(summary)[1]+1,]= c(chr_num,CRD_name,CRD_start,CRD_end,corr_coeff_mean_CRD_1,corr_coeff_mean_CRD_2,a@fisher1925$statistic,a@fisher1925$p.value)
  }
  summary
}

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94



for (chr_num in 1:22){
  
  # TODO: do we have to do it per chr?

  # Load files ---------------------------------------------

  
  CORR_neut = fread(file.path(path_peak, paste0(type, "_",  "neut", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE)
  CORR_mono = fread(file.path(path_peak, paste0(type, "_",  "mono", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE)
  CORR_tcell = fread(file.path(path_peak, paste0(type, "_",  "tcell", ".corr.chr", chr_num, ".txt")),
                     head = FALSE, stringsAsFactors = FALSE)
  
  MODS_neut = fread(file.path(path, paste0(type, "_", "neut", ".chr", chr_num, ".module.txt.gz")),
                    head = TRUE, stringsAsFactors = FALSE)
  MODS_mono = fread(file.path(path, paste0(type, "_", "mono", ".chr", chr_num, ".module.txt.gz")),
                    head = TRUE, stringsAsFactors = FALSE)
  MODS_tcell = fread(file.path(path, paste0(type, "_", "tcell", ".chr", chr_num, ".module.txt.gz")),
                     head = TRUE, stringsAsFactors = FALSE)
  
  # limit in cis the corr
  CORR_neut <- subset(CORR_neut, CORR_neut$V2 - CORR_neut$V1 <= 250)
  CORR_mono <- subset(CORR_mono, CORR_mono$V2 - CORR_mono$V1 <= 250)
  CORR_tcell <- subset(CORR_tcell, CORR_tcell$V2 - CORR_tcell$V1 <= 250)
  
  # only CRDs
  MODS_neut <- MODS_only_CRDs(MODS_neut)
  MODS_mono <- MODS_only_CRDs(MODS_mono)
  MODS_tcell <- MODS_only_CRDs(MODS_tcell)
  
  # compute mean correlation coeff per CRD   ---------------------------------------------
  # MODSs from cell type 1
  # CORR_2_c1 from cell type 1
  # CORR_2_c2 from cell type 2
  
  summary <- compute_corr_coeff_CRDs_and_significance(CORR_mono,CORR_neut, MODS_mono) 
  write.table(summary, file =file.path(outpath, paste0("corr_coeff_CRDs_chr",chr_num,"_","mono_vs_neut",".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  
  summary <- compute_corr_coeff_CRDs_and_significance(CORR_mono,CORR_tcell, MODS_mono) 
  write.table(summary, file =file.path(outpath, paste0("corr_coeff_CRDs_chr",chr_num,"_","mono_vs_tcell",".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  
  summary <- compute_corr_coeff_CRDs_and_significance(CORR_neut,CORR_tcell, MODS_neut) 
  write.table(summary, file =file.path(outpath, paste0("corr_coeff_CRDs_chr",chr_num,"_","neut_vs_tcell",".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  
  
  summary <- compute_corr_coeff_CRDs_and_significance(CORR_neut,CORR_mono, MODS_neut) 
  write.table(summary, file =file.path(outpath, paste0("corr_coeff_CRDs_chr",chr_num,"_","neut_vs_mono",".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  
  summary <- compute_corr_coeff_CRDs_and_significance(CORR_tcell,CORR_mono, MODS_tcell) 
  write.table(summary, file =file.path(outpath, paste0("corr_coeff_CRDs_chr",chr_num,"_","tcell_vs_mono",".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  
  summary <- compute_corr_coeff_CRDs_and_significance(CORR_tcell,CORR_neut, MODS_tcell) 
  write.table(summary, file =file.path(outpath, paste0("corr_coeff_CRDs_chr",chr_num,"_","tcell_vs_neut",".txt")), sep = "\t",
              row.names = FALSE, ,quote=FALSE)
  
}

 print("finised") 

COR_c1 <- CORR_neut
COR_c2 <- CORR_mono
MOD_cell_1 <- MODS_neut
crd_i <- 1
chr_num <- 1



