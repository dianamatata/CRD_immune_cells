# Clean environment ---------------------------------

rm(list=ls())
gc()


# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)

# Directories ---------------------------------

directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/conditional_merged/'
plot_directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/plots_pval_after_cond/'


# Main Loop ---------------------------------

file_list=list.files(path=directory,pattern=c('txt'))

for (file in file_list) {
  name=substr(file, 1, nchar(file)-7)
  df = as.data.frame(fread(paste0(directory,file),header=F))
  # Hist P val ---------------------------------
  
  qobj = tryCatch({ qvalue(df$V18)},
                error = function(error_condition) {
                  qvalue(df$V18, pi0=1)} )

#  qobj=qvalue(df$V18)
  pi0=qobj$pi0
  pi1=1-pi0

  print(paste0(name,' pi1 ', pi1))
  pdf(paste0(plot_directory,"histogram_p_",name,".pdf"))
  hist(df$V18, main=paste0(' pi0=',round(pi0,digits = 3),'  pi1=',round(pi1,digits = 3)), xlab='p values', breaks=20,cex.main=1.5,  cex.lab=1.2,cex.axis=1.2)
  dev.off()

}


