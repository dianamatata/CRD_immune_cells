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
outpath = "/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/PLOTS"


# Functions  ---------------------------------


MODS_per_section <- function(MODS, section) {
  MODS1 = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ] # otherwise too many triangles
  MODS2 <- MODS1 %>% filter(RIDX > section & LIDX < section+1000) # create a subset of MODS.
  MODS2
}

MODS_filtered <- function(MODS) {
  MODS1 = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ] # otherwise too many triangles
  MODS1
}



# Main ---------------------------------------------

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94
chr_num <- 5
section <- 1000



for (chr_num in 1:22){
  
  # Load files ---------------------------------------------
  
  MODS_neut = fread(file.path(path, paste0(type, "_", "neut", ".chr", chr_num, ".module.txt.gz")),
                    head = TRUE, stringsAsFactors = FALSE, select = c("LIDX","RIDX","UUID","N_REG","COUNT","MOD"))
  MODS_mono = fread(file.path(path, paste0(type, "_", "mono", ".chr", chr_num, ".module.txt.gz")),
                    head = TRUE, stringsAsFactors = FALSE, select = c("LIDX","RIDX","UUID","N_REG","COUNT","MOD"))
  MODS_tcell = fread(file.path(path, paste0(type, "_", "tcell", ".chr", chr_num, ".module.txt.gz")),
                     head = TRUE, stringsAsFactors = FALSE, select = c("LIDX","RIDX","UUID","N_REG","COUNT","MOD"))
  
  
  CORR_neut = fread(file.path(path_peak, paste0(type, "_",  "neut", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  CORR_mono = fread(file.path(path_peak, paste0(type, "_",  "mono", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  CORR_tcell = fread(file.path(path_peak, paste0(type, "_",  "tcell", ".corr.chr", chr_num, ".txt")),
                     head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  
  
  COR_mono_vs_neut = fread(file.path(path_corr, paste0('COR_mono_vs_neut_chr',chr_num,'_sorted.txt')),
                    head = FALSE, stringsAsFactors = FALSE)
  COR_mono_vs_tcell = fread(file.path(path_corr,  paste0('COR_mono_vs_tcell_chr',chr_num,'_sorted.txt')),
                           head = FALSE, stringsAsFactors = FALSE)
  COR_neut_vs_tcell = fread(file.path(path_corr,  paste0('COR_neut_vs_tcell_chr',chr_num,'_sorted.txt')),
                           head = FALSE, stringsAsFactors = FALSE)
  
  
  # limit in cis the corr  ---------------------------------------------
  
  CORR_neut <- subset(CORR_neut, CORR_neut$V2 - CORR_neut$V1 <= 250)
  CORR_mono <- subset(CORR_mono, CORR_mono$V2 - CORR_mono$V1 <= 250)
  CORR_tcell <- subset(CORR_tcell, CORR_tcell$V2 - CORR_tcell$V1 <= 250)
  
  MODS_neut <- MODS_filtered(MODS_neut)
  MODS_mono <- MODS_filtered(MODS_mono)
  MODS_tcell <- MODS_filtered(MODS_tcell)
  
  ## only significant ones
  
  pval <- 0.00005
  for (pval in c(0.00005, 0.0005, 0.005,0.05)){
    print(pval)
    
    COR_mono_vs_neut_F <-COR_mono_vs_neut %>% filter(V6 <= pval)
    COR_mono_vs_tcell_F <-COR_mono_vs_tcell %>% filter(V6 <= pval)
    COR_neut_vs_tcell_F <-COR_neut_vs_tcell %>% filter(V6 <= pval)
    
    # COR_mono_vs_neut_F$corrdiff=abs(COR_mono_vs_neut_F$V3)-abs(COR_mono_vs_neut_F$V4)
    # COR_mono_vs_tcell_F$corrdiff=abs(COR_mono_vs_tcell_F$V3)-abs(COR_mono_vs_tcell_F$V4)
    # COR_neut_vs_tcell_F$corrdiff=abs(COR_neut_vs_tcell_F$V3)-abs(COR_neut_vs_tcell_F$V4)
    COR_mono_vs_neut_F$corrdiff=abs(COR_mono_vs_neut_F$V3-COR_mono_vs_neut_F$V4)/2
    COR_mono_vs_tcell_F$corrdiff=abs(COR_mono_vs_tcell_F$V3-COR_mono_vs_tcell_F$V4)/2
    COR_neut_vs_tcell_F$corrdiff=abs(COR_neut_vs_tcell_F$V3-COR_neut_vs_tcell_F$V4)/2
    
    
    # Plot all chr  ---------------------------------------------
    
    # for one section, 250units in height, is 750 in print size
    # and width is 1500 units,  is 750*3 in print size
    rowplots <- 6
    ylimit_print <- 750*rowplots #*3
    xlimit_print <- max(CORR_mono$V2) #*3
    xlimit <- max(CORR_mono$V2)

    
    png(paste0(outpath,"/plot_chr_",chr_num,"_pval_",pval,"_.png"), width = xlimit_print, height = ylimit_print)
    par(mfrow=c(6,1),mar = c(3,3, 3, 3))
    
    ## MONO
    plot((CORR_mono$V2+CORR_mono$V1)/2, CORR_mono$V2-CORR_mono$V1, pch=18, col=rgb(0,0,1,abs(CORR_mono$V7)^2),  xlab=" ", ylab=" ",
         xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_mono)) {
      polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
    
    
    ## NEUT
    plot((CORR_neut$V2+CORR_neut$V1)/2, CORR_neut$V2-CORR_neut$V1, pch=18, col=rgb(0,0,1,abs(CORR_neut$V7)^2),  xlab=" ", ylab=" ",
         xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_neut)) {
      polygon(c(MODS_neut$LIDX[m], (MODS_neut$LIDX[m]+MODS_neut$RIDX[m])/2, MODS_neut$RIDX[m]), c(0,MODS_neut$RIDX[m]-MODS_neut$LIDX[m], 0), border="black", lwd=2)}
    
    ## TCELL
    plot((CORR_tcell$V2+CORR_tcell$V1)/2, CORR_tcell$V2-CORR_tcell$V1, pch=18, col=rgb(0,0,1,abs(CORR_tcell$V7)^2), xlab=" ", ylab=" ",
         xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_tcell)) {
      polygon(c(MODS_tcell$LIDX[m], (MODS_tcell$LIDX[m]+MODS_tcell$RIDX[m])/2, MODS_tcell$RIDX[m]), c(0,MODS_tcell$RIDX[m]-MODS_tcell$LIDX[m], 0), border="black", lwd=2)}
    
    
    plot((COR_mono_vs_neut_F$V2+COR_mono_vs_neut_F$V1)/2, COR_mono_vs_neut_F$V2-COR_mono_vs_neut_F$V1, pch=18, col=rgb(1,0,0,abs(COR_mono_vs_neut_F$corrdiff)^2),  xlab=" ", ylab=" ",
         xlim=c(1,xlimit), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_mono)) {
      polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
    for (m in 1:nrow(MODS_neut)) {
      polygon(c(MODS_neut$LIDX[m], (MODS_neut$LIDX[m]+MODS_neut$RIDX[m])/2, MODS_neut$RIDX[m]), c(0,MODS_neut$RIDX[m]-MODS_neut$LIDX[m], 0), border="grey", lwd=2)}
    
    plot((COR_mono_vs_tcell_F$V2+COR_mono_vs_tcell_F$V1)/2, COR_mono_vs_tcell_F$V2-COR_mono_vs_tcell_F$V1, pch=18, col=rgb(1,0,0,abs(COR_mono_vs_tcell_F$corrdiff)^2), xlab=" ", ylab=" ",
         xlim=c(1,xlimit), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_mono)) {
      polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
    for (m in 1:nrow(MODS_tcell)) {
      polygon(c(MODS_tcell$LIDX[m], (MODS_tcell$LIDX[m]+MODS_tcell$RIDX[m])/2, MODS_tcell$RIDX[m]), c(0,MODS_tcell$RIDX[m]-MODS_tcell$LIDX[m], 0), border="grey", lwd=2)}
    
    plot((COR_neut_vs_tcell_F$V2+COR_neut_vs_tcell_F$V1)/2, COR_neut_vs_tcell_F$V2-COR_neut_vs_tcell_F$V1, pch=18, col=rgb(1,0,0,abs(COR_neut_vs_tcell_F$corrdiff)^2), xlab=" ", ylab=" ",
         xlim=c(1,xlimit), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_neut)) {
      polygon(c(MODS_neut$LIDX[m], (MODS_neut$LIDX[m]+MODS_neut$RIDX[m])/2, MODS_neut$RIDX[m]), c(0,MODS_neut$RIDX[m]-MODS_neut$LIDX[m], 0), border="grey", lwd=2)}
    for (m in 1:nrow(MODS_tcell)) {
      polygon(c(MODS_tcell$LIDX[m], (MODS_tcell$LIDX[m]+MODS_tcell$RIDX[m])/2, MODS_tcell$RIDX[m]), c(0,MODS_tcell$RIDX[m]-MODS_tcell$LIDX[m], 0), border="grey", lwd=2)}
  }  
    dev.off()
    
    
}

#  
#  ## MONO
#  png(paste0(outpath,"/plot_chr_",chr_num,"_pval_",pval,"2.png"), xlimit_print, height = 750*3)
#  par(mfrow=c(1,1),mar = c(20,20, 20, 20))
# # plot(1, type="n", xlab="", ylab="",xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=20, cex.lab=20, cex.axis=10)
#  plot((CORR_mono$V2+CORR_mono$V1)/2, CORR_mono$V2-CORR_mono$V1, pch=18, col=rgb(0,0,1,abs(CORR_mono$V7)^2), xlab=" ", ylab=" ",
#       xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=20, cex.lab=20, cex.axis=10)
#  # for (m in 1:nrow(MODS_mono)) {
#  #   polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
#   dev.off()
#  
#  
# after compare_corr_coeff in which we compute
# Clean environment ------------------------------------

rm(list = ls())
gc()


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
outpath = "/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/PLOTS"


# Functions  ---------------------------------


MODS_per_section <- function(MODS, section) {
  MODS1 = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ] # otherwise too many triangles
  MODS2 <- MODS1 %>% filter(RIDX > section & LIDX < section+1000) # create a subset of MODS.
  MODS2
}

MODS_filtered <- function(MODS) {
  MODS1 = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ] # otherwise too many triangles
  MODS1
}



# Main ---------------------------------------------

type="hist"
Nneut <- 165
Nmono <- 160
Ntcell <- 94
chr_num <- 5
section <- 1000



for (chr_num in 1:22){
  
  # Load files ---------------------------------------------
  
  MODS_neut = fread(file.path(path, paste0(type, "_", "neut", ".chr", chr_num, ".module.txt.gz")),
                    head = TRUE, stringsAsFactors = FALSE, select = c("LIDX","RIDX","UUID","N_REG","COUNT","MOD"))
  MODS_mono = fread(file.path(path, paste0(type, "_", "mono", ".chr", chr_num, ".module.txt.gz")),
                    head = TRUE, stringsAsFactors = FALSE, select = c("LIDX","RIDX","UUID","N_REG","COUNT","MOD"))
  MODS_tcell = fread(file.path(path, paste0(type, "_", "tcell", ".chr", chr_num, ".module.txt.gz")),
                     head = TRUE, stringsAsFactors = FALSE, select = c("LIDX","RIDX","UUID","N_REG","COUNT","MOD"))
  
  
  CORR_neut = fread(file.path(path_peak, paste0(type, "_",  "neut", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  CORR_mono = fread(file.path(path_peak, paste0(type, "_",  "mono", ".corr.chr", chr_num, ".txt")),
                    head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  CORR_tcell = fread(file.path(path_peak, paste0(type, "_",  "tcell", ".corr.chr", chr_num, ".txt")),
                     head = FALSE, stringsAsFactors = FALSE,select = c(1,2,3,7))
  
  
  COR_mono_vs_neut = fread(file.path(path_corr, paste0('COR_mono_vs_neut_chr',chr_num,'_sorted.txt')),
                    head = FALSE, stringsAsFactors = FALSE)
  COR_mono_vs_tcell = fread(file.path(path_corr,  paste0('COR_mono_vs_tcell_chr',chr_num,'_sorted.txt')),
                           head = FALSE, stringsAsFactors = FALSE)
  COR_neut_vs_tcell = fread(file.path(path_corr,  paste0('COR_neut_vs_tcell_chr',chr_num,'_sorted.txt')),
                           head = FALSE, stringsAsFactors = FALSE)
  
  
  # limit in cis the corr  ---------------------------------------------
  
  CORR_neut <- subset(CORR_neut, CORR_neut$V2 - CORR_neut$V1 <= 250)
  CORR_mono <- subset(CORR_mono, CORR_mono$V2 - CORR_mono$V1 <= 250)
  CORR_tcell <- subset(CORR_tcell, CORR_tcell$V2 - CORR_tcell$V1 <= 250)
  
  MODS_neut <- MODS_filtered(MODS_neut)
  MODS_mono <- MODS_filtered(MODS_mono)
  MODS_tcell <- MODS_filtered(MODS_tcell)
  
  ## only significant ones
  
  pval <- 0.00005
  for (pval in c(0.00005, 0.0005, 0.005,0.05)){
    print(pval)
    
    COR_mono_vs_neut_F <-COR_mono_vs_neut %>% filter(V6 <= pval)
    COR_mono_vs_tcell_F <-COR_mono_vs_tcell %>% filter(V6 <= pval)
    COR_neut_vs_tcell_F <-COR_neut_vs_tcell %>% filter(V6 <= pval)
    
    # COR_mono_vs_neut_F$corrdiff=abs(COR_mono_vs_neut_F$V3)-abs(COR_mono_vs_neut_F$V4)
    # COR_mono_vs_tcell_F$corrdiff=abs(COR_mono_vs_tcell_F$V3)-abs(COR_mono_vs_tcell_F$V4)
    # COR_neut_vs_tcell_F$corrdiff=abs(COR_neut_vs_tcell_F$V3)-abs(COR_neut_vs_tcell_F$V4)
    COR_mono_vs_neut_F$corrdiff=abs(COR_mono_vs_neut_F$V3-COR_mono_vs_neut_F$V4)/2
    COR_mono_vs_tcell_F$corrdiff=abs(COR_mono_vs_tcell_F$V3-COR_mono_vs_tcell_F$V4)/2
    COR_neut_vs_tcell_F$corrdiff=abs(COR_neut_vs_tcell_F$V3-COR_neut_vs_tcell_F$V4)/2
    
    
    # Plot all chr  ---------------------------------------------
    
    # for one section, 250units in height, is 750 in print size
    # and width is 1500 units,  is 750*3 in print size
    rowplots <- 6
    ylimit_print <- 750*rowplots #*3
    xlimit_print <- max(CORR_mono$V2) #*3
    xlimit <- max(CORR_mono$V2)

    
    png(paste0(outpath,"/plot_chr_",chr_num,"_pval_",pval,"_.png"), width = xlimit_print, height = ylimit_print)
    par(mfrow=c(6,1),mar = c(3,3, 3, 3))
    
    ## MONO
    plot((CORR_mono$V2+CORR_mono$V1)/2, CORR_mono$V2-CORR_mono$V1, pch=18, col=rgb(0,0,1,abs(CORR_mono$V7)^2),  xlab=" ", ylab=" ",
         xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_mono)) {
      polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
    
    
    ## NEUT
    plot((CORR_neut$V2+CORR_neut$V1)/2, CORR_neut$V2-CORR_neut$V1, pch=18, col=rgb(0,0,1,abs(CORR_neut$V7)^2),  xlab=" ", ylab=" ",
         xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_neut)) {
      polygon(c(MODS_neut$LIDX[m], (MODS_neut$LIDX[m]+MODS_neut$RIDX[m])/2, MODS_neut$RIDX[m]), c(0,MODS_neut$RIDX[m]-MODS_neut$LIDX[m], 0), border="black", lwd=2)}
    
    ## TCELL
    plot((CORR_tcell$V2+CORR_tcell$V1)/2, CORR_tcell$V2-CORR_tcell$V1, pch=18, col=rgb(0,0,1,abs(CORR_tcell$V7)^2), xlab=" ", ylab=" ",
         xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_tcell)) {
      polygon(c(MODS_tcell$LIDX[m], (MODS_tcell$LIDX[m]+MODS_tcell$RIDX[m])/2, MODS_tcell$RIDX[m]), c(0,MODS_tcell$RIDX[m]-MODS_tcell$LIDX[m], 0), border="black", lwd=2)}
    
    
    
    plot((COR_mono_vs_neut_F$V2+COR_mono_vs_neut_F$V1)/2, COR_mono_vs_neut_F$V2-COR_mono_vs_neut_F$V1, pch=18, col=rgb(1,0,0,abs(COR_mono_vs_neut_F$corrdiff)^2),  xlab=" ", ylab=" ",
         xlim=c(1,xlimit), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_mono)) {
      polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
    for (m in 1:nrow(MODS_neut)) {
      polygon(c(MODS_neut$LIDX[m], (MODS_neut$LIDX[m]+MODS_neut$RIDX[m])/2, MODS_neut$RIDX[m]), c(0,MODS_neut$RIDX[m]-MODS_neut$LIDX[m], 0), border="grey", lwd=2)}
    
    plot((COR_mono_vs_tcell_F$V2+COR_mono_vs_tcell_F$V1)/2, COR_mono_vs_tcell_F$V2-COR_mono_vs_tcell_F$V1, pch=18, col=rgb(1,0,0,abs(COR_mono_vs_tcell_F$corrdiff)^2), xlab=" ", ylab=" ",
         xlim=c(1,xlimit), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_mono)) {
      polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
    for (m in 1:nrow(MODS_tcell)) {
      polygon(c(MODS_tcell$LIDX[m], (MODS_tcell$LIDX[m]+MODS_tcell$RIDX[m])/2, MODS_tcell$RIDX[m]), c(0,MODS_tcell$RIDX[m]-MODS_tcell$LIDX[m], 0), border="grey", lwd=2)}
    
    plot((COR_neut_vs_tcell_F$V2+COR_neut_vs_tcell_F$V1)/2, COR_neut_vs_tcell_F$V2-COR_neut_vs_tcell_F$V1, pch=18, col=rgb(1,0,0,abs(COR_neut_vs_tcell_F$corrdiff)^2), xlab=" ", ylab=" ",
         xlim=c(1,xlimit), ylim=c(0, 250),bty="n", cex.main=3, cex.lab=3, cex.axis=1)
    for (m in 1:nrow(MODS_neut)) {
      polygon(c(MODS_neut$LIDX[m], (MODS_neut$LIDX[m]+MODS_neut$RIDX[m])/2, MODS_neut$RIDX[m]), c(0,MODS_neut$RIDX[m]-MODS_neut$LIDX[m], 0), border="grey", lwd=2)}
    for (m in 1:nrow(MODS_tcell)) {
      polygon(c(MODS_tcell$LIDX[m], (MODS_tcell$LIDX[m]+MODS_tcell$RIDX[m])/2, MODS_tcell$RIDX[m]), c(0,MODS_tcell$RIDX[m]-MODS_tcell$LIDX[m], 0), border="grey", lwd=2)}
    
    dev.off()
    
    
  }
}

#  
#  ## MONO
#  png(paste0(outpath,"/plot_chr_",chr_num,"_pval_",pval,"2.png"), xlimit_print, height = 750*3)
#  par(mfrow=c(1,1),mar = c(20,20, 20, 20))
# # plot(1, type="n", xlab="", ylab="",xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=20, cex.lab=20, cex.axis=10)
#  plot((CORR_mono$V2+CORR_mono$V1)/2, CORR_mono$V2-CORR_mono$V1, pch=18, col=rgb(0,0,1,abs(CORR_mono$V7)^2), xlab=" ", ylab=" ",
#       xlim=c(1,xlimit ), ylim=c(0, 250),bty="n", cex.main=20, cex.lab=20, cex.axis=10)
#  # for (m in 1:nrow(MODS_mono)) {
#  #   polygon(c(MODS_mono$LIDX[m], (MODS_mono$LIDX[m]+MODS_mono$RIDX[m])/2, MODS_mono$RIDX[m]), c(0,MODS_mono$RIDX[m]-MODS_mono$LIDX[m], 0), border="black", lwd=2)}
#   dev.off()
#  
#  

