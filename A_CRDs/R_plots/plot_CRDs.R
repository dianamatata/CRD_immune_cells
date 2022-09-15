
# Clean environment ------------------------------------

rm(list=ls())
gc()

# Packages ---------------------------------------------

library(qvalue)
library(RColorBrewer)
library(data.table) # call fread
library(R.utils) # to read .gz
library("dplyr")

# path mac
path='/Users/dianaavalos/Programming/A_CRD_plots/1_CRD'
path_peak='/Users/dianaavalos/Programming/A_CRD_plots/8_PEAKS'

# path baobab to be updated

path='/home/users/a/avalosma/scratch/1_CRD'
path_peak='/home/users/a/avalosma/scratch/8_PEAKS'

type='methyl'
for (cell in c('neut','mono','tcell')){
  for (chr_num in 1:22){#22
    
    # get data
    file=paste0(type,"_",cell,".chr",chr_num,".module.txt.gz")
    MODS = fread(file.path(path,file), head=TRUE, stringsAsFactors=FALSE)
    
    file_corr=paste0(type,"_",cell,".corr.chr",chr_num,".txt")
    CORR2 = fread(file.path(path_peak,file_corr), head=FALSE, stringsAsFactors=FALSE)
    
    # limit in cis the corr
    CORR2$DIFF=CORR2$V2-CORR2$V1
    CORR2 <- subset(CORR2, DIFF <= 250)
    
    max_x=as.integer(max(CORR2$V2)/1000)*1000
    for (section in seq(1,max_x-1000,1000)) {
      print(paste0(type,"_",cell,"_chr",chr_num,' section ',section, ' ', section+1000))
      
      
      # plot computing
      MODSs = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ]
      MODS_FILTERED=MODS[MODS$MOD == 1 & MODS$N_REG > 1, ]
      
      #name_file
      filename=paste0(type,"_",cell,".chr",chr_num,"_",section,"_figCRD.png")
      outpath="/home/users/a/avalosma/scratch/A_CRD_plots/CRD_plots"
      
      # V6 became V2, V12: V8, #V7 corr value, V8 pval
      png(paste0(outpath,"/",filename), width=3000, height = 1000, res = 200, type="cairo") # with limits
      plot((CORR2$V2+CORR2$V1)/2, CORR2$V2-CORR2$V1, pch=18, col=rgb(0,0,1,abs(CORR2$V7)^2), xlab="corr", xlim=c(section, section +1000), ylim=c(0, 250)) #, bty="n", ylab="", xlab="Chromatin peak index on chromosome 14", yaxt="n")
      abline(h=0)
      for (m in 1:nrow(MODSs)) {
        polygon(c(MODSs$LIDX[m], (MODSs$LIDX[m]+MODSs$RIDX[m])/2, MODSs$RIDX[m]), c(0,MODSs$RIDX[m]-MODSs$LIDX[m], 0), border="black", lwd=2)
      }
      dev.off()
      
    }
  }
}







