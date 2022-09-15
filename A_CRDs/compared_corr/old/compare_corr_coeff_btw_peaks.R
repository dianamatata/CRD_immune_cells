
# Packages ---------------------------------------------

library(cocor)
library(dplyr)
library(data.table) # call fread
library(R.utils) # to read .gz

args = commandArgs(trailingOnly=TRUE)

# Paths ---------------------------------------------

# path mac
folder = '/Users/dianaavalos/Programming/A_CRD_plots/8_PEAKS'
# path baobab to be updated
folder = '/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/files'

# Function ---------------------------------------------

compared_corr<-function(C1,C2,C3)
{
  Nneut <- 165
  Nmono <- 160
  Ntcell <- 94
  
  a <- cocor.indep.groups(C1,C2, Nmono, Nneut, alternative = "two.sided",
                     test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                     data.name = NULL, var.labels = NULL, return.htest = FALSE)
  b <- cocor.indep.groups(C1,C3, Nmono, Ntcell, alternative = "two.sided",
                          test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                          data.name = NULL, var.labels = NULL, return.htest = FALSE)
  c <- cocor.indep.groups(C2,C3, Nneut, Ntcell, alternative = "two.sided",
                          test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                          data.name = NULL, var.labels = NULL, return.htest = FALSE)
  c(
    a@fisher1925$statistic,
    a@fisher1925$p.value,
    b@fisher1925$statistic,
    b@fisher1925$p.value,
    c@fisher1925$statistic,
    c@fisher1925$p.value
  )
}

# Main ---------------------------------------------
chr_num=args[1]
print(chr_num)
#chr_num <- 22
#for (chr_num in seq(1,22)){
print(chr_num)
CORR = fread(file.path(folder, paste0("mono_neut_tcell.chr", chr_num, ".txt")), header=FALSE)
f <- mapply(compared_corr,CORR$V5,CORR$V7,CORR$V9) # we need vectors
df <- data.frame(t(f))
colnames(df) <- c("mono_neut_rho","mono_neut_pval","mono_tcell_rho","mono_tcell_pval","neut_tcell_rho","neut_tcell_pval")
colnames(CORR) <- c("peak_i","peak_j","location_i","location_j","rho_mono","pval_mono","rho_neut","pval_neut","rho_tcell","pval_tcell")
COR <- cbind(CORR,df)
write.table(COR, file.path(folder,paste0("CORR_mono_neut_tcell_chr_",chr_num,".txt")), sep="\t", row.names=FALSE)
#}

print("finished")


