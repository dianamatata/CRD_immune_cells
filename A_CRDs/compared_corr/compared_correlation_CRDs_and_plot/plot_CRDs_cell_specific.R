# plot CRDs color cell specific



# Packages ---------------------------------------------

library(dplyr)
library(data.table) # call fread
library(R.utils) # to read .gz


# Paths ---------------------------------------------

# path mac
folder <- '/Users/dianaavalos/Programming/A_CRD_plots/13_CRD_corr'

# path baobab to be updated
folder <-  '/srv/beegfs/scratch/users/a/avalosma/13_CRD_corr'

chr=1
cell="mono"

for (chr in seq(1,22)){
  for( cell in c('mono','neut','tcell')){
    
    list_of_files <- list.files(path = folder, pattern = paste0("corr_coeff_CRDs_chr",chr,"_",cell))
    MODS <- fread(file.path(folder,list_of_files[1]))
    MODS2 <- fread(file.path(folder,list_of_files[2]))
    cell_top <- strsplit(strsplit(list_of_files[1],"[.]")[[1]][1],"[_]")[[1]][7]
    cell_bt <- strsplit(strsplit(list_of_files[2],"[.]")[[1]][1],"[_]")[[1]][7]
    
    # shade CRD= 1-pvalue so darker=more cell specific
    #png(paste0(path,"/plot_chr_",chr,"_",cell,"_top_",cell_top,".png"), width=3000, height = 1000, res = 200) #  type="cairo",with limits
    #plot(1, type="n", xlab="", ylab="", xlim=c(0, max(MODS$CRD_end)), ylim=c(-250, 250))
    
    rowplots <- 2
    ylimit_print <- 750*rowplots #*3
    xlimit_print <- max(MODS$CRD_end) #*3
    xlimit <- max(MODS$CRD_end)
    
    png(paste0(folder,"/plot_chr_",chr,"_",cell,"_top_",cell_top,".png"), width = xlimit_print, height = ylimit_print)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, xlimit), ylim=c(-250, 250),xaxt='n',yaxt='n')
    abline(h=0)
    axis(1, at=seq(0,max(MODS$CRD_end),1000),labels=seq(0,max(MODS$CRD_end),1000), col.axis="black", las=2)
    
    for (m in 1:nrow(MODS)) {
      polygon(c(MODS$CRD_start[m], (MODS$CRD_start[m]+MODS$CRD_end[m])/2, MODS$CRD_end[m]), c(0,MODS$CRD_end[m]-MODS$CRD_start[m], 0), col=rgb(0,0,1,1-MODS$pvalue[m]),border="black", lwd=1)
      polygon(c(MODS$CRD_start[m], (MODS$CRD_start[m]+MODS$CRD_end[m])/2, MODS$CRD_end[m]), c(0,-(MODS$CRD_end[m]-MODS$CRD_start[m]), 0), col=rgb(0,0,1,1-MODS2$pvalue[m]),border="black", lwd=1)
      
      # text(x=(MODS$CRD_start[m]+MODS$CRD_end[m])/2, y=MODS$CRD_end[m]-MODS$CRD_start[m]+10,labels=paste0(MODS$CRD_name[m],"_",m),col = "black",cex = 1)
    }
      dev.off()

  }
}


