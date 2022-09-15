
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
library(data.table)
library(GenomicRanges)
library(tidyverse)
library(stringr)

#############################################################################################
#
# FUNCTIONS
#
#############################################################################################


read.cormat <- function(filename,peakrange=250){
  cat(paste0("Importing ",filename," ...\n"))
  mydat.cor <- fread(filename)[,c(1,2,3,4,5,6,8)]
  colnames(mydat.cor) = c("ID1","ID2","pos1","pos2","PeakName1","PeakName2","pval")
  mydat.cor$peakrange=mydat.cor$ID2-mydat.cor$ID1
  mydat.cor$distance = mydat.cor$pos2-mydat.cor$pos1
  mydat.cor
}


compute_ratio_hic_mapdata <-function(CRDmindist,CRDmaxdist,cutoff=5,pval_cutoff=0.01){
  up = sum(abs(cordata$distance)>=CRDmindist & abs(cordata$distance)<CRDmaxdist & cordata_validated>cutoff & cordata$pval.adj>pval_cutoff,na.rm=T)
  down = sum(abs(cordata$distance)>=CRDmindist & abs(cordata$distance)<CRDmaxdist & cordata$pval.adj>pval_cutoff,na.rm=T)
  up/down*100
}

compute_ratio_hic_mapdata_signif <-function(CRDmindist,CRDmaxdist,cutoff=5,pval_cutoff=0.01){
  up = sum(abs(cordata$distance)>=CRDmindist & abs(cordata$distance)<CRDmaxdist & cordata_validated>cutoff & cordata$pval.adj<pval_cutoff,na.rm=T)
  down = sum(abs(cordata$distance)>=CRDmindist & abs(cordata$distance)<CRDmaxdist & cordata$pval.adj<pval_cutoff,na.rm=T)
  up/down*100
}

#############################################################################################
#
# MAIN
#
#############################################################################################
# directory='/Users/dianaavalos/Programming/A_CRD_plots/8_PEAKS'
# out_dir=directory
# PCHiC = fread('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/PCHiC_peak_matrix_cutoff5.tsv')


# for cluster
directory='/srv/beegfs/scratch/users/a/avalosma/8_PEAKS'
out_dir=directory
PCHiC = fread('/srv/beegfs/scratch/users/a/avalosma/0_CRD/PCHiC_peak_matrix_cutoff5.tsv')

buffersize = 500 #size of buffer around peak position midpoint
reImport = FALSE


colnames(PCHiC)[1] = "baitChr"
data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
data_type='hist'
cell_type='neut'
for(data_type in data_types){
  for(cell_type in cell_types){
    
    name_condition=paste0(data_type,'_',cell_type )
    cat(name_condition, '  \n')
    filename=paste0(out_dir,'/',name_condition,'_10R.txt')
    file.create(filename)
    
    #############################################################################################
    
    # import corr data
    all.files <- intersect(intersect(list.files(path=directory, pattern=data_type, full.names=TRUE, recursive=FALSE),
                                     list.files(path=directory, pattern=cell_type, full.names=TRUE, recursive=FALSE)),
                           list.files(path=directory, pattern="corr.chr", full.names=TRUE, recursive=FALSE))
    
    for (file in all.files){
      cat(file)
      
      l <- lapply(file, read.cormat)
      
      cordata <- rbindlist(l)
      
      cat("Generating adjusted pvalues...\n")
      cordata$pval.adj = p.adjust(cordata$pval,method='fdr')
      setkeyv(cordata,c("PeakName1","PeakName2"))
      
      # cordata$rowid = c(1:nrow(cordata))
      cordata$chr = str_match(cordata$PeakName1, "chr([0-9]*)")[,2]
      
      #############################################################################################
      
      # HiC
      baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
      oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))
      
      peak1bed <- GRanges(seqnames=cordata$chr,ranges=IRanges(start=cordata$pos1-buffersize, end=cordata$pos1+buffersize))
      peak2bed <- GRanges(seqnames=cordata$chr,ranges=IRanges(start=cordata$pos2-buffersize, end=cordata$pos2+buffersize))
      
      cat("Overlapping correlation and PCHiC data...\n")
      
      #fwd
      x = findOverlaps(baitbed,peak1bed)
      y = findOverlaps(oebed,peak2bed)
      tmp = rbind(as.data.frame(x),as.data.frame(y))
      validated.fwd = tmp[which(duplicated(tmp)),]
      
      #bwd
      x = findOverlaps(baitbed,peak2bed)
      y = findOverlaps(oebed,peak1bed)
      tmp = rbind(as.data.frame(x),as.data.frame(y))
      validated.bwd = tmp[which(duplicated(tmp)),]
      
      validated = unique(rbind(validated.fwd,validated.bwd))
      
      
      cat("Annotating correlation data...\n")
      
      cordata_validated = rep(0,nrow(cordata))
      
      for(i in 1:nrow(validated)){
        currenthic = validated$queryHits[i]
        currentmap = validated$subjectHits[i]
        currentHiCScore = PCHiC$Neu[currenthic]
        if(currentHiCScore>cordata_validated[currentmap]){
          cordata_validated[currentmap] = currentHiCScore
        }
      }
      
      cat("Annotating PCHiC data...\n")
      
      hic_validated = rep(NA,nrow(PCHiC))
      
      for(i in 1:nrow(validated)){
        currenthic = validated$queryHits[i]
        currentpval = cordata$pval.adj[validated$subjectHits[i]]
        if(currentpval<hic_validated[currenthic] | is.na(hic_validated[currenthic])){
          hic_validated[currenthic] = currentpval
        }
      }
      
      ### carefull cell type HiC
      currentHiCScore_all <- c(PCHiC$Neu[currenthic], PCHiC$Mon[currenthic], PCHiC$nCD4[currenthic])
      names(currentHiCScore_all) <- c("neut", "mono", "tcell")
      currentHiCScore=currentHiCScore_all[[cell_type]]
      
      # myhist_bg = hist(PCHiC$Neu[hic_validated<1],breaks = c(0,5,10,15,20,2000),plot=F)
      # myhist_signif = hist(PCHiC$Neu[hic_validated<0.01],breaks = c(0,5,10,15,20,2000),plot=F)
      # 
      myhist_bg = hist(currentHiCScore[hic_validated<1],breaks = c(0,5,10,15,20,2000),plot=F)
      myhist_signif = hist(currentHiCScore[hic_validated<0.01],breaks = c(0,5,10,15,20,2000),plot=F)
      write(file,file=filename,append=TRUE)
      write(myhist_bg,file=filename,append=TRUE)
      write(myhist_signif,file=filename,append=TRUE)
      
      
      
      ############################################### FIG
      
      pdf(paste0(out_dir,'/',name_condition,"_HiC_validation.pdf"),5,5)
      toplot = data.frame(counts = myhist_signif$counts/myhist_bg$counts*100,Number = c("0-5","5-10","10-15","15-20",">20"))
      toplot$Number = factor(toplot$Number,levels = c("0-5","5-10","10-15","15-20",">20"))
      g <- ggplot(toplot, aes(x = Number, y = counts))+ ggtitle("HiC contacts with CRD associations") +
        geom_bar(stat = "identity",fill="#E69F00") +
        labs(x = "PC HiC Score",y = "Fraction of correlated peak pairs (%)") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
      print(g)
      dev.off()
      
      dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06)
      write(toplot,file=filename,append=TRUE)
      write(toplot$Number,file=filename,append=TRUE)
      
      
      mat_hic = matrix(0,nrow=(length(dist_bins)-1),ncol=1)
      mat_hic_signif = matrix(0,nrow=(length(dist_bins)-1),ncol=1)
      
      for(i in 1:(length(dist_bins)-1)){
        mat_hic[i,1] = compute_ratio_hic_mapdata(dist_bins[i],dist_bins[i+1])
        mat_hic_signif[i,1] = compute_ratio_hic_mapdata_signif(dist_bins[i],dist_bins[i+1])
      }
      
      ############################################### FIG
      
      pdf(paste0(out_dir,'/',name_condition,"_HiC_support_gene_CRD_associations.pdf"))
      toplot = data.frame(correlated = mat_hic_signif[,1], uncorrelated = mat_hic[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
      toplot$dist <- factor(toplot$dist ,levels = c("0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
      dfm <- melt(toplot,id.vars = 3)
      
      ggplot(dfm, aes(x = dist, y = value))+ ggtitle("HiC support for peak associations") +
        geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +
        labs(x = "Distance",y = "Fraction with PCHiC support (%)") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))
      
      dev.off()
      
      save.image(paste0(out_dir,'/',name_condition,"_Correlation_vs_PCHiC.RData"))
      
      stop()
      
      ############################################### FIG
      
      
      pdf(paste0(out_dir,'/',name_condition,"_scatterplot.pdf"))
      plot(-log10(cordata$pval[which(abs(cordata$distance>1e3)==1 & cordata_validated>0)]),cordata_validated[which(abs(cordata$distance>1e3)==1 & cordata_validated>0)])
      plot(-log10(cordata$pval[which(abs(cordata$distance>1e4)==1 & cordata_validated>0)]),cordata_validated[which(abs(cordata$distance>1e4)==1 & cordata_validated>0)])
      plot(-log10(cordata$pval[which(abs(cordata$distance>1e5)==1 & cordata_validated>0)]),cordata_validated[which(abs(cordata$distance>1e5)==1 & cordata_validated>0)])
      dev.off()
      
      cor.test(-log10(cordata$pval[which((abs(cordata$distance>0))==1 & cordata_validated>0)]),cordata_validated[which((abs(cordata$distance>0))==1 & cordata_validated>0)],method='s')
      
    }        
  }
}









