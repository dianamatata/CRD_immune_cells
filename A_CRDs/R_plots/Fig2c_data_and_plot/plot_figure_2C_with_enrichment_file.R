rm(list = ls())

# plot figure 2C of CRD paper 
library(data.table)    
library(dplyr)
library(tidyverse)
library(tidyselect)
library(ggpubr)
library(grid)
library(gridExtra)

# PLOT GARFIELD RESULTS
summary <- read.csv(file=paste0("/Users/dianaavalos/Desktop/reviews_avalos/CRD_immune_cells/A_CRDs/R_plots/Fig2c_data_and_plot/CRDQTL_outfiles/enrichment_hCRD_QTLs_plot2c.txt"), sep="\t", header=TRUE,skip = 0)
thresh=1e-05
sumCRDQTLs <- summary[summary$PThresh==thresh,]
sumCRDQTLs <- sumCRDQTLs[sumCRDQTLs$GWAS %in% c("BC","EC", "LC","MC","NT", "PC","RBC","RC" ,"WBC","Height"),]
sumCRDQTLs$GWAS <- factor(sumCRDQTLs$GWAS,levels = c("BC","EC", "LC","MC","NT", "PC","RBC","RC" ,"WBC","Height"))

# PLOT
listPlotCRDQTLs1 = list()
i=1
Annot="MON hCRD"
for (Annot in unique(sumCRDQTLs$Annotation)[1:2]){
  print(Annot)
  dff <- sumCRDQTLs[sumCRDQTLs$Annot==Annot,]
  listPlotCRDQTLs1[[i]] <-
    ggplot(dff) +
    geom_errorbar( aes(x=GWAS, ymin=SEmin, ymax=SEmax), width=0.2, colour="black", alpha=1, size=0.5) +
    geom_bar( aes(x=GWAS, y=OR, fill=logpval), stat="identity", alpha=1, width = 0.8) + theme_minimal() +
    scale_fill_gradient(low = "#0066FF", high = "#FF3300",limits = c(0,ceiling(max(summary$logpval))))+
    theme(axis.text.x = element_text(size=0,color="white")) + labs(x = " ", y = paste0(Annot)) + scale_y_continuous(breaks=c(1,2,3,4)) +ylim(0,4)
  i=i+1
}
Annot=unique(sumCRDQTLs$Annotation)[3]
dff <- sumCRDQTLs[sumCRDQTLs$Annot==Annot,]
listPlotCRDQTLs1[[i]] <-
  ggplot(dff) +
  geom_errorbar( aes(x=GWAS, ymin=SEmin, ymax=SEmax), width=0.2, colour="black", alpha=1, size=0.5) +
  geom_bar( aes(x=GWAS, y=OR, fill=logpval), stat="identity", alpha=1, width = 0.8) + theme_minimal() +
  scale_fill_gradient(low = "#0066FF", high = "#FF3300",limits = c(0,ceiling(max(summary$logpval))))+
  theme(axis.text.x = element_text(size=10,color="black")) + labs(x = " ", y = paste0(Annot)) + scale_y_continuous(breaks=c(1,2,3,4)) +ylim(0,4)

png(paste0("figure_CRD_blood_GWAS_",thresh,".png"), 2100, 2100, res=300)
grid.arrange(grobs = listPlotCRDQTLs1, nrow = length(unique(sumCRDQTLs$Annotation))+1)
dev.off()


