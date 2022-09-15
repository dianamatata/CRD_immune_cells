library(data.table)
library(ggplot2)

# ran in computer directly
directory='/home/users/a/avalosma/scratch/4_CRD_residualized'
outdirectory='/home/users/a/avalosma/scratch/A_CRD_plots/QTLs_PCs_plots'
#directory='/Users/dianaavalos/Programming/A_CRD_plots/QTLs_PCs_plots'
summary_QTLs=fread(paste0(directory,'/summary_significant_cisQTLs.txt'))
for (cell in unique(summary_QTLs$cell_type)){
  for (type in unique(summary_QTLs$folder)){
    print(paste0(cell,' ', type))
    subset_cell=dplyr::filter(summary_QTLs, grepl(cell, cell_type))
    subset_cell_type=dplyr::filter(subset_cell, grepl(type, folder))
    title=paste0(outdirectory,'/PCs_QTLs_discovered_',cell,'_',type,'.png')
    print(length(subset_cell_type$PC))
    png(title, width = 2000, height = 2000, res=300)
    plot(subset_cell_type$PC,subset_cell_type$significant_QTLs,type = "b",
         main = paste0(type, '-QTLs discovered in ',cell),
         xlab = "PCs",
         ylab = paste0(type,"-QTLs discovered") )
    dev.off()
  }
}

