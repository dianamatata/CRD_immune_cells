### old!!!!! and not working

# Clean environment
rm(list=ls())
gc()


library(qvalue)
library(ggplot2)
library(gplots)
library(data.frame)
library(GenomicRanges)
library(RColorBrewer)
library("circlize")
library(igraph)

#### which ones for tcells? CD4 not CD8!!!
PCHiC = fread('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/PCHiC_peak_matrix_cutoff5.tsv')
PCHiC_Neu=PCHiC$Neu
PCHiC_Mon=PCHiC$Mon
PCHiC_Tcl=PCHiC$tCD4

#######  get PCHiC data

colnames(PCHiC)[1] = "baitChr"
interchromosomal = which(PCHiC$baitChr!=PCHiC$oeChr)
PCHiC = PCHiC[interchromosomal,]
baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))
CRD1bed <- GRanges(seqnames=trans_data$chr1,ranges=IRanges(start=trans_data$start1, end=trans_data$end1))
CRD2bed <- GRanges(seqnames=trans_data$chr2,ranges=IRanges(start=trans_data$start2, end=trans_data$end2))

#fwd
x = findOverlaps(baitbed,CRD1bed)
y = findOverlaps(oebed,CRD2bed)
tmp = rbind(as.trans_data.frame(x),as.trans_data.frame(y))
validated.fwd = tmp[which(duplicated(tmp)),]

#bwd
x = findOverlaps(baitbed,CRD2bed)
y = findOverlaps(oebed,CRD1bed)
tmp = rbind(as.trans_data.frame(x),as.trans_data.frame(y))
validated.bwd = tmp[which(duplicated(tmp)),]

validated = unique(rbind(validated.fwd,validated.bwd))

hic_validated = rep(1,nrow(PCHiC))

for(i in 1:nrow(validated)){
  currenthic = validated$queryHits[i]
  currentqval = trans_data$qval[validated$subjectHits[i]]
  if(currentqval<hic_validated[currenthic]){
    hic_validated[currenthic] = currentqval
  }
}

####* loop around files and label cell types

path = "/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant_bis/"
file.names <- dir(path, pattern =".txt")

for(i in 1:length(file.names)){
  cat(file.names[i], '  ')
  filename=file.names[i]
  name= substr(file.names[i],1,nchar(file.names[i])-21)
  
  if (str_detect(name, 'tcell')){
    PCHiC_cell=PCHiC_Tcl
  }
  if (str_detect(name, 'mono')){
    PCHiC_cell=PCHiC_Mon
  }
  if (str_detect(name, 'neut')){
    PCHiC_cell=PCHiC_Neu
  }
  
  
  #######  WRITE SIGNIFICANT HITS
  
  trans_data = as.data.frame(data.table::fread(paste0(path,'/',filename), head=FALSE, stringsAsFactors=FALSE))
  # weird number of columns
  colnames(trans_data) = c("idx1","chr1","start1","end1","id1","idx2","chr2","start2","end2","id2","corr","pval")
  Q = qvalue(trans_data$pval) # keep trans_data when qvalue < 0.01 test
  trans_data$qval = Q$qvalues
  trans_data_q001 = trans_data[Q$qvalue < 0.01, ]
  
  write.table(trans_data_q001, paste0(name,"_test.significant1FDR.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  
  ###### 5.5 5e ?
  
  myhist_bg = hist(PCHiC_cell,breaks = c(0,10,20,50,100,2000),plot=F)
  myhist_signif = hist(PCHiC_cell[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)
  
  pdf(paste0(name,"_HiC_validation.pdf"),5,5)
  toplot = trans_data.frame(counts = myhist_signif$counts/myhist_bg$counts*100,Number = c("0-10","11-20","21-50","51-100",">100"))
  toplot$Number = factor(toplot$Number,levels = c("0-10","11-20","21-50","51-100",">100"))
  g <- ggplot(toplot, aes(x = Number, y = counts))+ ggtitle("HiC contacts with CRD associations") +
    geom_bar(stat = "identity",fill="#E69F00") +
    labs(x = "PC HiC Score",y = "Fraction (%)") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
  print(g)
  dev.off()
  
  
  trans_data_q001$chr1 = paste0("chr",trans_data_q001$chr1)
  trans_data_q001$chr2 = paste0("chr",trans_data_q001$chr2)
  
  
  #PLOT1: HISTOGRAM + PI1
  pdf(paste0(name,"_Histogram.pdf"),10,5)
  par(mfrow=c(1, 2))
  hist(trans_data$pval, xlab="Nominal P-values", main="", breaks=100)
  abline(h=Q$pi0 * nrow(trans_data) / 100, col="red")
  legend("topright", legend=c(paste("pi1=", signif(100-Q$pi0*100, 3), "%", sep=""), paste("#hits=", nrow(trans_data_q001), " (1% FDR)", sep="")), bty="n")
  hist(trans_data_q001$corr, breaks=100, main="", xlab="Correlation coefficient")
  hist(Q$qvalue,xlab="Q-values", main="", breaks=100)
  dev.off()
  
  #PLOT2: EXAMPLE HEATMAP
  #TD = as.trans_data.frame(trans_data.table::fread(paste("zcat ~/VitalIT/SGX/V2/trans_data_correlation/trans/plot/LCL_ALL.",C1,".",C2,".plot.txt.gz", sep="")), head=FALSE)
  #M1 = as.trans_data.frame(trans_data.table::fread(paste("zcat ~/VitalIT/SGX/V2/trans_data_correlation/cis/chromatin/LCL_ALL.chr",C1,".subset.txt.gz", sep="")), head=FALSE)
  #M2 = as.trans_data.frame(trans_data.table::fread(paste("zcat ~/VitalIT/SGX/V2/trans_data_correlation/cis/chromatin/LCL_ALL.chr",C2,".subset.txt.gz", sep="")), head=FALSE)
  
  #C1=10; C2=16; C1c=c(4000, 4500); C2c=c(3250, 3750);
  #TDs = TD[TD$V1 > C1c[1] & TD$V1 < C1c[2] & TD$V3 > C2c[1] & TD$V3 < C2c[2], ]
  #M1s = M1[M1$V1 > C1c[1] & M1$V5 < C1c[2], ]
  #M2s = M2[M2$V1 > C2c[1] & M2$V5 < C2c[2], ]
  
  #pdf("figureB.4.2.pdf", 8,8)
  #plot(TDs$V1-C1c[1], TDs$V3-C2c[1], xlim=c(-125, 500), ylim=c(-125, 500), xaxt="n", yaxt="n", xlab=paste("Chromosome", C1), ylab=paste("Chromosome", C2), col=rgb(0,0,1, TDs$V5^2), pch=15, main="", cex=0.5)
  #abline(h=0); abline(v=0);
  #points((M1s$V1 + M1s$V5)/2-C1c[1], (M1s$V1-M1s$V5)*0.5, col=rgb(0,0,1, M1s$V10^2), pch=15, cex=0.5)
  #points((M2s$V1-M2s$V5)*0.5, (M2s$V1 + M2s$V5)/2-C2c[1], col=rgb(0,0,1, M2s$V10^2), pch=15, cex=0.5)
  #axis(1, at=seq(0, 500, 100), labels=seq(C1c[1], C1c[2], 100))
  #axis(2, at=seq(0, 500, 100), labels=seq(C2c[1], C2c[2], 100))
  #dev.off()
  
  #PLOT3: CIRCLE PLOT FOR BEST 1000 LINKS
  COL = brewer.pal(9,"Set1")
  trans_data_q001s = trans_data_q001[order(trans_data_q001$pval),]
  trans_data_q001s = trans_data_q001s[1:1000, ]
  
  pdf("Circle_plot.pdf", 8, 8)
  par(mar=c(1, 1, 1, 1))
  bed1 = trans_data.frame(chr=trans_data_q001s$chr1, from=trans_data_q001s$start1, to=trans_data_q001s$end1)
  bed2 = trans_data.frame(chr=trans_data_q001s$chr2, from=trans_data_q001s$start2, to=trans_data_q001s$end2)
  dir = (trans_data_q001s$corr > 0)
  circos.clear()
  circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 22:1))
  circos.genomicLink(bed1, bed2, col=rgb(0,0,1,0.05))
  dev.off()
  
  #PLOT4: CONNECTIVITY
  pdf(paste0(name,"_Connectivity.pdf"),6,6)
  hist(table(c(trans_data_q001$id1, trans_data_q001$id2)), breaks=40, xlab="Number of connected modules (module degree)", main="")
  conhist = hist(table(c(trans_data_q001$id1, trans_data_q001$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  frac = conhist$counts/sum(conhist$counts)*100
  barplot(frac,names = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  toplot = trans_data.frame(counts = frac,Number = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  toplot$Number = factor(toplot$Number,levels = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  ggplot(toplot, aes(x = Number, y = counts))+ ggtitle("Connectivity of CRD trans associations") +
    geom_bar(stat = "identity",fill="#E69F00") +
    labs(x = "Connectivity",y = "Fraction (%)") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
  
  dev.off()
  
  #PLOT4: CONNECTIVITY 3 
  pdf(paste0(name,"_Connectivity2.pdf"),6,6)
  hist(table(c(trans_data_q001.NEU$id1, trans_data_q001.NEU$id2)), breaks=40, xlab="Number of connected modules (module degree)", main="")
  conhistNEU = hist(table(c(trans_data_q001.NEU$id1, trans_data_q001.NEU$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  conhistMON = hist(table(c(trans_data_q001.MON$id1, trans_data_q001.MON$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  conhistTCL = hist(table(c(trans_data_q001.TCL$id1, trans_data_q001.TCL$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  
  fracNEU = conhistNEU$counts/sum(conhistNEU$counts)*100
  fracMON = conhistMON$counts/sum(conhistMON$counts)*100
  fracTCL = conhistTCL$counts/sum(conhistTCL$counts)*100
  
  barplot(frac,names = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  toplot = trans_data.frame(counts = frac,Number = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  toplot$Number = factor(toplot$Number,levels = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  ggplot(toplot, aes(x = Number, y = counts))+ ggtitle("Connectivity of CRD trans associations") +
    geom_bar(stat = "identity",fill="#E69F00") +
    labs(x = "Connectivity",y = "Fraction (%)") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
  
  dev.off()
  
  
  #PLOT5: NETWORK EXAMPLE OF THE 1000 BEST LINKS
  library('igraph')
  library(RColorBrewer)
  COL1 = brewer.pal(8,"Set1")
  COL2 = brewer.pal(8,"Set2")
  COL3 = brewer.pal(8,"Set3")
  COL=c(COL1, COL2, COL3)
  REF = trans_data.frame(i=1:22, s=paste("", 1:22, sep=""))
  # trans_data with qvalue < 0.01, ordered by pvalue, take first 1000 nodes
  trans_data_q001s = trans_data_q001[order(trans_data_q001$pval),]
  trans_data_q001s = trans_data_q001s[1:1000, ]
  
  LINKS = trans_data.frame(from=trans_data_q001s$id1, to=trans_data_q001s$id2, weigth=ifelse(trans_data_q001s$corr>0, 1, 2), stringsAsFactors=FALSE)
  NODES = trans_data.frame(id=names(table(c(LINKS$from, LINKS$to))), chr=matrix(unlist(strsplit(names(table(c(LINKS$from, LINKS$to))), split="_")), ncol=3, byrow=TRUE)[, 1], stringsAsFactors=FALSE)
  trans_datanet = graph_from_trans_data_frame(d=LINKS, vertices=NODES, directed=FALSE)
  deg <- degree(trans_datanet, mode="all")
  
  # added by me
  write.table(LINKS, paste0(name,"_LINKS.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  V(trans_datanet)$size <- log2(deg)+1
  V(trans_datanet)$frame.color <- "white"
  V(trans_datanet)$color <- COL[match(V(trans_datanet)$chr, REF$s)]
  V(trans_datanet)$border <- "n"
  V(trans_datanet)$label <- ""
  E(trans_datanet)$arrow.mode <- 0
  E(trans_datanet)$color <- ifelse(E(trans_datanet)$weigth == 1, rgb(0,0,1, 0.4), rgb(1,0,0, 0.4))
  
  pdf(paste0(name,"_igraph_network.pdf"), 8, 8)
  par(mar=c(1,1,1,1))
  plot(trans_datanet)
  legend("topleft", fill=COL[1:22], legend=REF$s, bg="white", title="Nodes:", ncol=11, cex=0.8)
  legend("bottomleft", fill=c(rgb(0,0,1, 0.4), rgb(1,0,0, 0.4)), legend=c("Positively correlated", "Negatively correlated"), bg="white", title="Links:")
  dev.off()
  
  ###### PLOT6: CALL TRHs using greedy algorithm
  
  LINKS.TRH = trans_data.frame(from=trans_data_q001$id1, to=trans_data_q001$id2, weigth=ifelse(trans_data_q001$corr>0, 1, 2), stringsAsFactors=FALSE)
  NODES.TRH = trans_data.frame(id=names(table(c(LINKS.TRH$from, LINKS.TRH$to))), chr=matrix(unlist(strsplit(names(table(c(LINKS.TRH$from, LINKS.TRH$to))), split="_")), ncol=3, byrow=TRUE)[, 1], stringsAsFactors=FALSE)
  trans_datanet.TRH = graph_from_trans_data_frame(d=LINKS.TRH, vertices=NODES.TRH, directed=FALSE)
  communities  = fastgreedy.community(trans_datanet.TRH)
  N_TRH = max(communities$membership)
  N_TRH_TOSHOW = 50
  pdf(paste0(name,"_TRH_group_sizes.pdf"))
  df = trans_data.frame(ID=1:N_TRH,SIZE=rle(sort(communities$membership))$length,stringsAsFactors=F)
  df <-df[order(-df$SIZE),]
  df <- trans_data.frame(head(df$SIZE,50),seq(1,50))
  names(df)[1]='SIZE'
  names(df)[2]='ID'
  
  p2<- ggplot(df, aes(x = ID, y = SIZE)) + geom_bar(stat = "identity", fill="steelblue", width=0.7) + labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH), x ="TRH index", y = "TRH Size") 
  p2 + theme(axis.text=element_text(size=16), axis.title=element_text(size=16) ) + theme_minimal() + scale_y_continuous(trans = 'log10')
  
  #head(df,10)
  
  #### g version
  pdf(paste0(name,"_TRH_group_sizes.pdf"))
  p<- ggplot(df,aes(x = ID, y = SIZE))
  p+labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH), x ="TRH index", y = "TRH Size") +
    geom_bar(stat="identity", fill="steelblue", width=0.7)+
    theme_minimal() + xlim(0, N_TRH_TOSHOW) + scale_y_continuous(trans = 'log10') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))
  print(p)
  dev.off()
  
  
  
  #PLOT7: CALL chromosome pair frequencies
  
  trans_data_q001ig = trans_data[trans_data$pval<1e-06,]
  chromFreq = matrix(0,22,22)
  for(i in 1:21){
      for(j in (i+1):22){
          frac_signif = sum(trans_data_q001ig$chr1 == i & trans_data_q001ig$chr2 == j)/sum(trans_data$chr1 == i & trans_data$chr2 == j)
          chromFreq[i,j] = frac_signif
          chromFreq[j,i] = frac_signif
      }
  }
  
  # chromFreq = sweep(chromFreq,MARGIN=1,FUN="/",STATS=rowSums(chromFreq))
  coul <- colorRampPalette(brewer.pal(9, "Blues"))(15)
  
  pdf(paste0(name,"_Frequency_by_chrom_pairs.pdf"))
  heatmap.2(chromFreq,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col=coul)
  dev.off()
  
  save(PCHiC,trans_data,trans_data_q001,trans_data_q001s,hic_validated,Q,file=paste0(name,"_TRANS_analysis.rda"))

}