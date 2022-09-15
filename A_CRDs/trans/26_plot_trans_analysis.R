# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(qvalue)
library(ggplot2)
library(gplots)
library(data.frame)
library(GenomicRanges)
library(RColorBrewer)
library("circlize")
library(igraph)


#############################################################################################
#
# INTEGRATION OF HIC WITH TRANS DATA by with compute_hic_validated
#
#############################################################################################

compute_hic_validated <- function(PCHiC, TRH_signif){
  
  colnames(PCHiC)[1] = "baitChr"
  interchromosomal = which(PCHiC$baitChr!=PCHiC$oeChr)
  PCHiC = PCHiC[interchromosomal,]
  
  baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
  oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))
  
  CRD1bed <- GRanges(seqnames=TRH_signif$chr1,ranges=IRanges(start=TRH_signif$start1, end=TRH_signif$end1))
  CRD2bed <- GRanges(seqnames=TRH_signif$chr2,ranges=IRanges(start=TRH_signif$start2, end=TRH_signif$end2))
  
  #fwd
  x = findOverlaps(baitbed,CRD1bed)
  y = findOverlaps(oebed,CRD2bed)
  tmp = rbind(as.data.frame(x),as.data.frame(y))
  validated.fwd = tmp[which(duplicated(tmp)),]
  
  #bwd
  x = findOverlaps(baitbed,CRD2bed)
  y = findOverlaps(oebed,CRD1bed)
  tmp = rbind(as.data.frame(x),as.data.frame(y))
  validated.bwd = tmp[which(duplicated(tmp)),]
  
  validated = unique(rbind(validated.fwd,validated.bwd))
  
  hic_validated = rep(1,nrow(PCHiC))
  
  for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentqval = TRH_signif$qval[validated$subjectHits[i]]
    if(currentqval<hic_validated[currenthic]){
      hic_validated[currenthic] = currentqval
    }
  }
  hic_validated
}


#############################################################################################
#
# Fig 5C FUNCTION
#
#############################################################################################

plot_Connectivity_ALL_CellTypes <- function(transCRD_NEU,transCRD_MON,transCRD_TCL, name){
  
  pdf(paste0(out_dir,"5c_Connectivity_ALL_CellTypes_",name,".pdf"), 6, 6)
  hist.NEU = hist(table(c(transCRD_NEU$id1, transCRD_NEU$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  hist.MON = hist(table(c(transCRD_MON$id1, transCRD_MON$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  hist.TCL = hist(table(c(transCRD_TCL$id1, transCRD_TCL$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  
  frac.NEU = hist.NEU$counts/sum(hist.NEU$counts)*100
  frac.MON = hist.MON$counts/sum(hist.MON$counts)*100
  frac.TCL = hist.TCL$counts/sum(hist.TCL$counts)*100
  
  toplot = data.frame(NEU = frac.NEU,MON = frac.MON,TCL = frac.TCL,Number = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  toplot.melted = reshape2::melt(toplot,id.vars = "Number", measure.vars = c("NEU","MON","TCL"))
  colnames(toplot.melted)[2] = "CellType"
  toplot.melted$Number = factor(toplot$Number,levels = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  
  ### now plot Fig 5C
  g <- ggplot(toplot.melted, aes(x = Number, y = value,fill=CellType))+ ggtitle("") +
    geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07" )) +
    labs(x = "Connectivity",y = "Fraction of CRDs (%)") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
  
  print(g)
  dev.off()
}

## if max bigger than 500 in breaks ? try to debug
# 
# hist.NEU = hist(table(c(transCRD_NEU$id1, transCRD_NEU$id2)), breaks=c(0,5,10,20,50,100,200,max(max(table(c(transCRD_NEU$id1, transCRD_NEU$id2))))),plot=F)
# hist.MON = hist(table(c(transCRD_MON$id1, transCRD_MON$id2)), breaks=c(0,5,10,20,50,100,200,max(max(table(c(transCRD_NEU$id1, transCRD_NEU$id2))))),plot=F)
# hist.TCL = hist(table(c(transCRD_TCL$id1, transCRD_TCL$id2)), breaks=c(0,5,10,20,50,100,200,max(max(table(c(transCRD_NEU$id1, transCRD_NEU$id2))))),plot=F)
#   

#############################################################################################
#
# Fig 5E FUNCTION
#
#############################################################################################

plot_transCRD_HiC <- function(toplot, name){
  
  pdf(paste0(out_dir,"5e_transCRD_HiC_ALL_CellTypes_",name,".pdf"), 5, 7)
  
  toplot.melted = reshape2::melt(toplot,id.vars = "Number", measure.vars = c("NEU","MON","TCL"))
  colnames(toplot.melted)[2] = "CellType"
  toplot.melted$Number = factor(toplot$Number,levels = c("0-10","11-20","21-50","51-100",">100"))
  g <- ggplot(toplot.melted, aes(x = Number, y = value,fill=CellType))+ ggtitle("") +
    geom_bar(position="dodge", stat="identity") + labs(x = "PC HiC Score",y = "Fraction of associations (%)") +  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07" )) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
  print(g)
  
  dev.off()
}


#############################################################################################
#
# COMPUTE transCRD association with PCHiC per cell type
#
#############################################################################################

compute_transCRDassociation_PCHiC <- function(PCHiC, transCRD_NEU,transCRD_MON,transCRD_TCL){
  
  hic_validated_NEU=compute_hic_validated(PCHiC, transCRD_NEU) # TRH_signif=transCRD_NEU
  hic_validated_MON=compute_hic_validated(PCHiC, transCRD_MON)
  hic_validated_TCL=compute_hic_validated(PCHiC, transCRD_TCL)
  
  myhist_bg.NEU = hist(PCHiC$Neu,breaks = c(0,10,20,50,100,2000),plot=F)
  myhist_signif.NEU = hist(PCHiC$Neu[hic_validated_NEU<0.05],breaks = c(0,10,20,50,100,2000),plot=F)
  
  myhist_bg.MON = hist(PCHiC$Mon,breaks = c(0,10,20,50,100,2000),plot=F)
  myhist_signif.MON = hist(PCHiC$Mon[hic_validated_MON<0.05],breaks = c(0,10,20,50,100,2000),plot=F)
  
  myhist_bg.TCL = hist(PCHiC$nCD4,breaks = c(0,10,20,50,100,2000),plot=F)
  myhist_signif.TCL = hist(PCHiC$nCD4[hic_validated_TCL<0.05],breaks = c(0,10,20,50,100,2000),plot=F)
  
  toplot = data.frame(Number = c("0-10","11-20","21-50","51-100",">100"),
                      NEU = myhist_signif.NEU$counts/myhist_bg.NEU$counts*100, 
                      MON = myhist_signif.MON$counts/myhist_bg.MON$counts*100,
                      TCL = myhist_signif.TCL$counts/myhist_bg.TCL$counts*100)
  toplot
}
  # toplot5
  # myhist_bg.MON$counts is similar, myhist_signif.MON$counts is diff. because hic_validated_MON and hic_validated diff
  # should have 603  49  39  22  18, 
  # I have 526  34  28   9   5
  


# topoplot 5
# Number      NEU      MON      TCL
# 1   0-10 14.45989 13.32008 13.29442
# 2  11-20 12.88136 14.45428 14.69534
# 3  21-50 20.44199 15.60000 18.39623
# 4 51-100 34.00000 29.72973 34.66667
# 5   >100 43.24324 37.50000 30.43478

# toplot
# Number      NEU      MON      TCL
# 1   0-10 14.01070 11.61917 12.55945
# 2  11-20 11.52542 10.02950 13.26165
# 3  21-50 19.33702 11.20000 16.50943
# 4 51-100 34.00000 12.16216 34.66667
# 5   >100 43.24324 10.41667 30.43478

#############################################################################################
#
# DIRECTORIES AND FILES
#
#############################################################################################

directory='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant'
directory2='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant_bis'

out_dir='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:other_plots/'

files <- list.files(path=directory, pattern="0.0*.txt", full.names=TRUE, recursive=FALSE)

# PCHiC$Neu , PCHiC$Mon, PCHiC$tCD4 #### which ones for tcells? CD4 not CD8!!!
PCHiC = fread('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/PCHiC_peak_matrix_cutoff5.tsv')


#############################################################################################
#
# MAIN
#
#############################################################################################

Paper_colors=c("#00AFBB", "#E7B800", "#FC4E07" )
FDR='0.05'

## plot 5c
for(data_type in c('hist','methyl')){
  for(module in c('mean','loom')){

      name=paste0(data_type,'_',module,'_',FDR)
      cat (name, '  ')
      
      files <- intersect(list.files(path=directory, pattern=data_type, full.names=TRUE, recursive=FALSE),
                         intersect(list.files(path=directory, pattern=module, full.names=TRUE, recursive=FALSE),
                                   list.files(path=directory, pattern=FDR, full.names=TRUE, recursive=FALSE)))
      
      transCRD_NEU = as.data.frame(data.table::fread(files[grepl("neu", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
      transCRD_MON = as.data.frame(data.table::fread(files[grepl("mon", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
      transCRD_TCL = as.data.frame(data.table::fread(files[grepl("tc", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
      #colnames(transCRD_NEU) = colnames(transCRD_MON) = colnames(transCRD_TCL) = c("idx1","chr1","start1","end1","id1","idx2","chr2","start2","end2","id2","corr","pval","qval","midplace","midplace2")
      
      a=try(plot_Connectivity_ALL_CellTypes(transCRD_NEU,transCRD_MON,transCRD_TCL, name))
      if(class(a) == "try-error"){print('error')}
       #plot_Connectivity_ALL_CellTypes_max_table(transCRD_NEU,transCRD_MON,transCRD_TCL, name)}
    
  }
}


data_type='hist'
module='mean'
FDR='0.01'

## plot 5e
for(data_type in c('hist','methyl')){
  for(module in c('mean','loom')){

    name=paste0(data_type,'_',module,'_',FDR)
    cat (name, '  ')
    
    files <- intersect(list.files(path=directory2, pattern=data_type, full.names=TRUE, recursive=FALSE),
                       intersect(list.files(path=directory2, pattern=module, full.names=TRUE, recursive=FALSE),
                                 list.files(path=directory2, pattern=FDR, full.names=TRUE, recursive=FALSE)))

    transCRD_NEU = as.data.frame(data.table::fread(files[grepl("neu", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
    transCRD_MON = as.data.frame(data.table::fread(files[grepl("mon", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
    transCRD_TCL = as.data.frame(data.table::fread(files[grepl("tce", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))

    # hic_validated=compute_hic_validated(PCHiC, TRH_signif) inside compute_transCRDassociation_PCHiC
    toplot2 <-try(compute_transCRDassociation_PCHiC(PCHiC, transCRD_NEU,transCRD_MON,transCRD_TCL))
    print(toplot2)
    if(class(toplot) != "try-error")
    {plot_transCRD_HiC(toplot2, name)}
    if(class(toplot) == "try-error")
    {print('err')}
  }
}



