# Overlap across cell types for chromatin peaks belonging to a CRD in a given cell type.


# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES and PATHS
#
#############################################################################################
library(UpSetR)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggmosaic)
Paper=c("#00AFBB", "#E7B800", "#FC4E07","#972D15")
color_palette= Paper

# basically grep gene column then comupte interest
directory='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_peaks/peaks/'
out_directory='Programming/A_CRD_plots/1_CRD/plots/'


#############################################################################################
#
# MAIN
#
#############################################################################################

data_types = list('hist','methyl')

data_type='hist'
condition='mean'

for(data_type in data_types){
  
  name=paste0(data_type,'_',condition)
  cat(name)
  
  NEU=fread(paste0(directory,data_type,'_','neut','_',condition,'.ALLchr.peaks.txt'), header = FALSE)$V1
  MON=fread(paste0(directory,data_type,'_','mono','_',condition,'.ALLchr.peaks.txt'),header = FALSE)$V1
  TCL=fread(paste0(directory,data_type,'_','tcell','_',condition,'.ALLchr.peaks.txt'),header = FALSE)$V1

  
  N_M = length(intersect(NEU,MON))
  N_T = length(intersect(NEU,TCL))
  M_T = length(intersect(TCL,MON))
  N_M_T = length(intersect(NEU,intersect(TCL,MON)))
  
  NEU1 = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
  MON1 = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
  TCL1 = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T)
  
  toplot2 = data.frame(NEU=c(length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+N_T-N_M_T+0,N_M_T),
                       MON=c(length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+0+M_T-N_M_T,N_M_T),
                       TCL=c(length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T),0+N_T-N_M_T+M_T-N_M_T,N_M_T),Group =c("1cell","2cells","3cells"))
  
  
  toplot2.molten = melt(toplot2,id.vars="Group")
  colnames(toplot2.molten)[2] = "CellType"
  toplot2.molten$Group <- factor(toplot2.molten$Group ,levels = c("1cell","2cells","3cells"))
  
  ############################################################################################# PLOT1
  expressionInput <- c(NEU = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), MON = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), 
                       TCL = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T), 
                       `NEU&MON` = N_M-N_M_T, `NEU&TCL` = N_T-N_M_T,`MON&TCL` = M_T-N_M_T, `NEU&MON&TCL` = N_M_T)
  
  pdf(paste0(out_directory,name,"_Overlap_chromatin_peaks.pdf"),paper='a4r')
  #upset(fromExpression(expressionInput), order.by = "degree",text.scale=1.8)
  upset(fromExpression(expressionInput), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")
  dev.off()
}
