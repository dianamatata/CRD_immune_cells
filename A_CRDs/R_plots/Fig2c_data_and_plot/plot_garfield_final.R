rm(list = ls())

library(data.table)    
library(dplyr)
library(tidyverse)
library(tidyselect)
library(ggpubr)
library(grid)
library(gridExtra)

# PLOT GARFIELD RESULTS

# CRDQTLs

dir="/Users/dianaavalos/Desktop/reviews_avalos/enrichment_QTLs/CRDQTL_outfiles"
for(pattern in c("_5.out")){
  print(pattern)
  files=list.files(dir,pattern)
  
  summary <- data.frame(matrix(ncol = 15, nrow = 0))
  for (f in files){
    dff <- fread(file.path(dir,f), select = c(1:14))
    GWAS <- str_split(f, '[.]')[[1]][1]
    dff$GWAS <- GWAS
    summary <- rbind(summary,dff)
  }
  
  # rename diseases and annotation
  summary$Annotation <- str_remove(summary$Annotation, "cisCRDQTLs_")
  summary$Annotation <- str_remove(summary$Annotation, "_mean.txt")
  summary$Annotation <- str_remove(summary$Annotation, ".txt")
  summary$GWAS <- str_remove(summary$GWAS, "BT_")
  summary$GWAS <- str_remove(summary$GWAS, "_10")
  summary$GWAS <- str_remove(summary$GWAS, "_5")
  summary$GWAS <- str_remove(summary$GWAS, "_5")
  
  summary[summary$Annotation=="hist_neut",]$Annotation="NEU hCRD"
  summary[summary$Annotation=="hist_mono",]$Annotation="MON hCRD"
  summary[summary$Annotation=="hist_tcell",]$Annotation="TCL hCRD"
  summary[summary$Annotation=="methyl_neut",]$Annotation="NEU mCRD"
  summary[summary$Annotation=="methyl_mono",]$Annotation="MON mCRD"
  summary[summary$Annotation=="methyl_tcell",]$Annotation="TCL mCRD"
  summary[summary$GWAS=="GIANT_HEIGHT",]$GWAS="Height"
  
  # values
  summary$logpval=-log10(summary$Pvalue)
  summary <- data.frame(summary)
  summary <- summary[summary$Annotation %in% unique(summary$Annotation)[1:6],]
  summary$SEmin <- summary$OR-summary$SE
  summary$SEmax <- summary$OR+summary$SE
  summary[summary$SEmin<=0.1,]$SEmin=0.1
  summary[summary$SEmax>4,]$SEmax=4
  
  # sumCRDQTLs[sumCRDQTLs$OR>10,]$OR=10
  # sumCRDQTLs[sumCRDQTLs$SE>1,]$SE=1
  thresh=1e-05
  print(thresh)
  sumCRDQTLs <- summary[summary$PThresh==thresh,]
  sumCRDQTLs$GWAS <- factor(sumCRDQTLs$GWAS,                                   
                            levels =  c("BC","EC", "LC","MC","NT", "PC","RBC","RC" ,"WBC", "RA" ,"CEL","MS","IBD","UC","CD" ,"DT1","DT2","Height"))
  
  # filter only BT
  # sumCRDQTLs <- sumCRDQTLs[sumCRDQTLs$GWAS %in% c("BC","EC", "LC","MC","NT", "PC","RBC","RC" ,"WBC","Height"),]
  # sumCRDQTLs$GWAS <- factor(sumCRDQTLs$GWAS,levels = c("BC","EC", "LC","MC","NT", "PC","RBC","RC" ,"WBC","Height"))
  # 
  # PLOT
  listPlotCRDQTLs1 = list()
  i=1
  Annot="MON hCRD"
  for (Annot in unique(sumCRDQTLs$Annotation)[1:5]){
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
  Annot=unique(sumCRDQTLs$Annotation)[6]
  dff <- sumCRDQTLs[sumCRDQTLs$Annot==Annot,]
  listPlotCRDQTLs1[[i]] <-
    ggplot(dff) +
    geom_errorbar( aes(x=GWAS, ymin=SEmin, ymax=SEmax), width=0.2, colour="black", alpha=1, size=0.5) +
    geom_bar( aes(x=GWAS, y=OR, fill=logpval), stat="identity", alpha=1, width = 0.8) + theme_minimal() +
    scale_fill_gradient(low = "#0066FF", high = "#FF3300",limits = c(0,ceiling(max(summary$logpval))))+
    theme(axis.text.x = element_text(size=10,color="black")) + labs(x = " ", y = paste0(Annot)) + scale_y_continuous(breaks=c(1,2,3,4)) +ylim(0,4)
  
  png(paste0("/Users/dianaavalos/Desktop/figure_CRD_blood_GWAS_ALL_",pattern,".png"), 2100, 2100, res=300)
  grid.arrange(grobs = listPlotCRDQTLs1, nrow = length(unique(sumCRDQTLs$Annotation))+1)
  dev.off()
  
}


