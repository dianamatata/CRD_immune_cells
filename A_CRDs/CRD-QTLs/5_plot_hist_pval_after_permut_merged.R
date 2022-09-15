# Clean environment ---------------------------------

rm(list=ls())
gc()


# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)

# Directories ---------------------------------

directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/merged_1000/'
plot_directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/plots_pval/'


# Main Loop ---------------------------------

file_list=list.files(path=directory,pattern=c('txt.gz'))

for (file in file_list) {
  name=substr(file, 1, nchar(file)-12)
  df = as.data.frame(fread(paste0(directory,file),header=F))
  # Hist P val ---------------------------------
  
  qobj=qvalue(df$V20)
  pi0=qobj$pi0
  pi1=1-pi0
  pdf(paste0(plot_directory,"histogram_p_",name,".pdf"))
  hist(df$V20, main=paste0(' pi0=',round(pi0,digits = 3),'  pi1=',round(pi1,digits = 3)), xlab='p values', breaks=20,cex.main=1.5,  cex.lab=1.2,cex.axis=1.2)
  dev.off()

}


