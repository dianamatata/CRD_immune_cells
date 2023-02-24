# Paths ---------------------------------------------
library(png)

wdir="/Users/dianaavalos/switchdrive/plot_CRD_genes_main_fig2"
setwd(wdir)
path=paste0(wdir,'/CRD')
path_peak=paste0(wdir,'/CORR')
path_associations=paste0(wdir,"/CRD_GENE")
outpath=paste0(wdir,"/OUT")
path_RNA=paste0(wdir,"/RNA")
rna_file <- c('EGAD00001002675_chr5.txt', 'EGAD00001002674_chr5.txt', 'EGAD00001002671_chr5.txt')
names(rna_file) <- c("neut", "mono", "tcell")
section=1000


  # for (chr_num in 1:22){#22
type="hist"
chr_num=5
  # LOAD  CRDs and peak correlation ---------------------------------------------
MODS_N = fread(file.path(path, paste0(type, "_", 'neut', ".chr", chr_num, ".module.txt.gz")),head = TRUE,stringsAsFactors = FALSE)
MODS_M = fread(file.path(path, paste0(type, "_", 'mono', ".chr", chr_num, ".module.txt.gz")),head = TRUE,stringsAsFactors = FALSE)
MODS_T = fread(file.path(path, paste0(type, "_", 'tcell', ".chr", chr_num, ".module.txt.gz")),head = TRUE,stringsAsFactors = FALSE)

CORR2_N = fread(file.path(path_peak, file_corr = paste0(type, "_", 'neut', ".corr.chr", chr_num, ".txt")),head = FALSE,stringsAsFactors = FALSE) %>% filter(V1 >=section -200 & V2<=section+1000+250)
CORR2_M = fread(file.path(path_peak, file_corr = paste0(type, "_", 'mono', ".corr.chr", chr_num, ".txt")),head = FALSE,stringsAsFactors = FALSE) %>% filter(V1 >=section -200 & V2<=section+1000+250)
CORR2_T = fread(file.path(path_peak, file_corr = paste0(type, "_", 'tcell', ".corr.chr", chr_num, ".txt")),head = FALSE,stringsAsFactors = FALSE) %>% filter(V1 >=section -200 & V2<=section+1000+250)

CORR2_M <- subset(CORR2_M, CORR2_M$V2-CORR2_M$V1 <= 200)
CORR2_N <- subset(CORR2_N, CORR2_N$V2-CORR2_N$V1 <= 200)
CORR2_T <- subset(CORR2_T, CORR2_T$V2-CORR2_T$V1 <= 200)

# middle_of_peak <- function (LSF) { 
#   LSF$pos1=(LSF$V3+LSF$V4)/2
#   LSF$pos2=(LSF$V8+LSF$V9)/2
#   LSF}
# i already have the middle
head(CORR2_N)


# LOAD  genes correlation ---------------------------------------------

# GENES_N=read.table(paste0("/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots","/txt/ALL_genes_100filter_","_",5,"hist","_",'neut',".chr",chr_num,"_",section,"2.txt"), sep = "\t", header=T)
# GENES_M=read.table(paste0("/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots","/txt/ALL_genes_100filter_","_",5,"hist","_",'mono',".chr",chr_num,"_",section,"2.txt"), sep = "\t", header=T)
# GENES_T=read.table(paste0("/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots","/txt/ALL_genes_100filter_","_",5,"hist","_",'tcell',".chr",chr_num,"_",section,"2.txt"), sep = "\t", header=T)
# GENES_N[GENES_N$associated==1,]$hgnc_symbol

CRD_count=4
for (CRD_count in c(1,3,4,5)) {
  print(CRD_count)
  
  GENES_N=read.table(paste0(outpath,"/txt/ALL_genes_100filter_","_",CRD_count,"hist","_","neut",".chr",chr_num,"_",section,"2.txt"), sep = "\t", header=T)
  GENES_M=read.table(paste0(outpath,"/txt/ALL_genes_100filter_","_",CRD_count,"hist","_","mono",".chr",chr_num,"_",section,"2.txt"), sep = "\t", header=T)
  GENES_T=read.table(paste0(outpath,"/txt/ALL_genes_100filter_","_",CRD_count,"hist","_","tcell",".chr",chr_num,"_",section,"2.txt"), sep = "\t", header=T)
  
  GENES_ALL= rbind(GENES_N,GENES_M,GENES_T)  %>% select(c(hgnc_symbol,start))
  GENES_ALL=unique(GENES_ALL)
  
  # modify MODS  ---------------------------------------------
  # before MODS =5
  
  MODSs_N <- MODS_N[MODS_N$N_REG > 1 & MODS_N$MOD == 1 & MODS_N$COUNT > CRD_count, ]  %>% filter(RIDX > section & LIDX < section+1000) # create a subset of MODS in the window 1000-2000
  MODSs_M <- MODS_M[MODS_M$N_REG > 1 & MODS_M$MOD == 1 & MODS_M$COUNT > CRD_count, ]  %>% filter(RIDX > section & LIDX < section+1000) # create a subset of MODS in the window 1000-2000
  MODSs_T <- MODS_T[MODS_T$N_REG > 1 & MODS_T$MOD == 1 & MODS_T$COUNT > CRD_count, ]  %>% filter(RIDX > section & LIDX < section+1000) # create a subset of MODS in the window 1000-2000
  
  GENES_M_notexpressed=GENES_ALL[!GENES_ALL$hgnc_symbol %in% GENES_M$hgnc_symbol,]
  GENES_N_notexpressed=GENES_ALL[!GENES_ALL$hgnc_symbol %in% GENES_N$hgnc_symbol,]
  GENES_T_notexpressed=GENES_ALL[!GENES_ALL$hgnc_symbol %in% GENES_T$hgnc_symbol,]
  
  
  i_min1=900;i_max1=i_min1+250;i_min2=2000;i_max2=2250; 
  p_min1 = round(CORR2_M$V3[CORR2_M$V1 == i_min1][1]) # min((CORR2_M$V3))
  p_max1 = round(CORR2_M$V4[CORR2_M$V2 == i_max1][1]) # 
  p_min2 = round(CORR2_M$V3[CORR2_M$V1 == i_min2][1])
  p_max2 = round(CORR2_M$V4[CORR2_M$V2 == i_max2][1])
  idx2pos1 <- function (i_curr) { unit_idx1 = (p_max2 - p_min1) / (i_max2 - i_min1); return (p_min1 + (i_curr - i_min1) * unit_idx1); }
  
  # PLOT ---------------------------------------------
  
  # cut before  49739081 which is EMB, so maybe 48 000 000 ?

    
  pdf(paste0("/Users/dianaavalos/Desktop/figure5_corrected_CRD_",CRD_count,".pdf"), 15, 10)
  
  # png(paste0("/Users/dianaavalos/Desktop/figure5_corrected_CRD_",CRD_count,".png"), 15, 10)
  # plot.new()
  
  
  Y_M=0
  Y_N=-400
  Y_T=-800
  par(fig=c(0.05, 1.00, 0.05, 1.00), new=FALSE) # c(x1, x2, y1, y2)
  plot(0,0, type="n", xlim=c(p_max1, p_min2), ylim=c(-980, 380), xaxt="n", yaxt="n", xlab="Genomic location on chromosome 5 (Mbp)", ylab="", main="Example of gene-CRD associations")
  box(lwd=1) 
  axis(1, at=seq(3e6, 15e6, 1e6), labels=seq(3, 15, 1))
  axis(2, at=c(125,-40,-80,-150), labels=c("Squared\ncorrelation\n","CRDs","TSSs","Genes"), las=2, tick=TRUE, lty=1, cex.axis=0.8)
  axis(2, at=c(125+Y_N,-40+Y_N,-80+Y_N,-150+Y_N), labels=c("Squared\ncorrelation\n","CRDs","TSSs","Genes"), las=2, tick=TRUE, lty=1, cex.axis=0.8)
  axis(2, at=c(125+Y_T,-40+Y_T,-80+Y_T,-150+Y_T), labels=c("Squared\ncorrelation\n","CRDs","TSSs","Genes"), las=2, tick=TRUE, lty=1, cex.axis=0.8)
  axis(1, at=c(40*1000000,45*1000000,50*1000000,55*1000000,60*1000000), labels=c(40,45,50,55,60), las=2, tick=TRUE, lty=1, cex.axis=0.8)
  abline(h=c(0,0+Y_T,0+Y_N), lwd=2)
  abline(h=c(-40-40+Y_N,-40+Y_T), col="lightgrey")
  abline(h=c(-80,-80+Y_N,-80+Y_T), col="lightgrey")
  
  
  
  for (l in 1:nrow(GENES_ALL)) {
    segments(GENES_ALL$start[l], 500, GENES_ALL$start[l],200, col="lightgrey",lty=1) # 2 dashed
    segments(GENES_ALL$start[l], -80, GENES_ALL$start[l], -200, col="lightgrey",lty=1)
    segments(GENES_ALL$start[l], Y_N-80, GENES_ALL$start[l],Y_N-200, col="lightgrey",lty=1)
    segments(GENES_ALL$start[l], Y_T-80, GENES_ALL$start[l],Y_T-200, col="lightgrey",lty=1)} 
  # add gene expression ticks
  log_gene_exp=c(-2,0,2)
  ticks=Y_N-80-(log_gene_exp+3)*20
  log_gene_exp_labels=c(0.01,1,100)
  axis(4, at=Y_M-80-(log_gene_exp+3)*20, labels=log_gene_exp_labels, las=2, tick=TRUE, lty=1, cex.axis=0.7)
  axis(4, at=Y_T-80-(log_gene_exp+3)*20, labels=log_gene_exp_labels, las=2, tick=TRUE, lty=1, cex.axis=0.7)
  axis(4, at=Y_N-80-(log_gene_exp+3)*20, labels=log_gene_exp_labels, las=2, tick=TRUE, lty=1, cex.axis=0.7)
  # axis(4, at=Y_M-80, labels="gene_exp", las=2, tick=TRUE, lty=1, cex.axis=0.7)
  
  # add gene infos
  GENES_ALL=GENES_ALL[order(GENES_ALL$start),]
  text(GENES_ALL$start[1]+100000, 380, GENES_ALL$hgnc_symbol[1],cex = 0.6, pos = 2, col = "grey30",srt = 90)
  for (l in 2:nrow(GENES_ALL)) {
    posy=400 
    if(l%%2 ==0){posy=285}
    text(GENES_ALL$start[l]+100000, posy, GENES_ALL$hgnc_symbol[l],cex = 0.6, pos = 2, col = "grey30",srt = 90) }# 220000 spacing}
  
  for (l in 1:nrow(GENES_M)) {
    if (GENES_M$associated[l]==1){
      segments(GENES_M$start[l], -80, GENES_M$start[l], -80-(GENES_M$logmean[l]+3)*20, col="#E7B800",lty=1,lwd=2)
      points(GENES_M$start[l], -80-(GENES_M$logmean[l]+3)*20, col="#E7B800", pch=20, cex=1)}
    if(GENES_M$associated[l]==0){
      segments(GENES_M$start[l], -80, GENES_M$start[l], -80-(GENES_M$logmean[l]+3)*20, col="azure4",lty=1,lwd=2)
      points(GENES_M$start[l], -80-(GENES_M$logmean[l]+3)*20, col="azure4", pch=20, cex=01) }}
  for (l in 1:nrow(GENES_N)) {
    if (GENES_N$associated[l]==1){segments(GENES_N$start[l], Y_N-80, GENES_N$start[l], Y_N-80-(GENES_N$logmean[l]+3)*20, col="#E7B800",lty=1,lwd=2)
      points(GENES_N$start[l], Y_N-80-(GENES_N$logmean[l]+3)*20, col="#E7B800", pch=20, cex=1)}
    if(GENES_N$associated[l]==0){
      segments(GENES_N$start[l], Y_N-80, GENES_N$start[l], Y_N-80-(GENES_N$logmean[l]+3)*20, col="azure4",lty=1,lwd=2)
      points(GENES_N$start[l], Y_N-80-(GENES_N$logmean[l]+3)*20, col="azure4", pch=20, cex=01) }}
  for (l in 1:nrow(GENES_T)) {
    if (GENES_T$associated[l]==1){segments(GENES_T$start[l], Y_T-80, GENES_T$start[l], Y_T-80-(GENES_T$logmean[l]+3)*20, col="#E7B800",lty=1,lwd=2)
      points(GENES_T$start[l], Y_T-80-(GENES_T$logmean[l]+3)*20, col="#E7B800", pch=20, cex=1)}
    if(GENES_T$associated[l]==0){
      segments(GENES_T$start[l], Y_T-80, GENES_T$start[l], Y_T-80-(GENES_T$logmean[l]+3)*20, col="azure4",lty=1,lwd=2)
      points(GENES_T$start[l], Y_T-80-(GENES_T$logmean[l]+3)*20, col="azure4", pch=20, cex=01) }}
  
  # add genes not in GENES_ALL
  for (l in 1:nrow(GENES_M_notexpressed)) {  points(GENES_M_notexpressed$start[l], -80, col="azure4", pch=20, cex=1)}
  for (l in 1:nrow(GENES_N_notexpressed)) {  points(GENES_N_notexpressed$start[l], Y_N-80, col="azure4", pch=20, cex=1)}
  for (l in 1:nrow(GENES_T_notexpressed)) {  points(GENES_T_notexpressed$start[l], Y_T-80, col="azure4", pch=20, cex=1)}
  
  
  # for GENEM signif or not # (GENES_N$logmean+3)*20
  
  
  # CRD MONO
  points(idx2pos1((CORR2_M$V1+CORR2_M$V2)/2), CORR2_M$V2-CORR2_M$V1, col=rgb(0,0,1,abs(CORR2_M$V7)^2), pch=18, cex=0.6)
  for (m in 1:nrow(MODSs_M)) {
    polygon(c(idx2pos1(MODSs_M$LIDX[m]), idx2pos1((MODSs_M$LIDX[m]+MODSs_M$RIDX[m])/2), idx2pos1(MODSs_M$RIDX[m])), c(0,MODSs_M$RIDX[m]-MODSs_M$LIDX[m], 0), border="black", lwd=2)}
  for (m in 1:nrow(MODSs_M)) {
    polygon(c(idx2pos1(MODSs_M$LIDX[m]), idx2pos1(MODSs_M$RIDX[m]), MODSs_M$END[m], MODSs_M$START[m]), c(0,0,-40,-40), border=NA, col=rgb(0,0,1, alpha=0.3))
    arrows(MODSs_M$START[m], -40, MODSs_M$END[m], -40, length=0.025, angle=90, code=3, col="blue")}
  for (l in 1:nrow(GENES_M)) {
    if(GENES_M$associated[l]==1){
      segments(GENES_M$start[l], -80, (GENES_M$V11[l]+GENES_M$V10[l])/2, -40, col="black") # V11 V10 start end crd
      points((GENES_M$V11[l]+GENES_M$V10[l])/2, -40, col="blue", pch=20, cex=0.5)}
    }
  for (g in 1:nrow(GENES_M)) {
    X = c(ifelse(GENES_M$V5[g] == "+", GENES_M$start[g]-2e4, GENES_M$stop[g]+2e4), ifelse(GENES_M$V5[g] == "+", GENES_M$start[g]-2e4, GENES_M$stop[g]+2e4), ifelse(GENES_M$V5[g] == "+", GENES_M$start[g]+8e4, GENES_M$stop[g]-8e4))
    polygon(X, c(-90, -70, -80), col="red", border=NA)}
  
  # CRD NEUT
  
  points(idx2pos1((CORR2_N$V1+CORR2_N$V2)/2), (CORR2_N$V2-CORR2_N$V1)+Y_N, col=rgb(0,0,1,abs(CORR2_N$V7)^2), pch=18, cex=0.6)
  for (m in 1:nrow(MODSs_N)) {
    polygon(c(idx2pos1(MODSs_N$LIDX[m]), idx2pos1((MODSs_N$LIDX[m]+MODSs_N$RIDX[m])/2), idx2pos1(MODSs_N$RIDX[m])), c(0+Y_N,(MODSs_N$RIDX[m]-MODSs_N$LIDX[m])+Y_N, +Y_N), border="black", lwd=2)}
  for (m in 1:nrow(MODSs_N)) {
    polygon(c(idx2pos1(MODSs_N$LIDX[m]), idx2pos1(MODSs_N$RIDX[m]), MODSs_N$END[m], MODSs_N$START[m]), c(0+Y_N,0+Y_N,-40+Y_N,-40+Y_N), border=NA, col=rgb(0,0,1, alpha=0.3))
    arrows(MODSs_N$START[m], -40+Y_N, MODSs_N$END[m], -40+Y_N, length=0.025, angle=90, code=3, col="blue")}
  for (l in 1:nrow(GENES_N)) {if(GENES_N$associated[l]==1){
    segments(GENES_N$start[l], -80+Y_N, (GENES_N$V11[l]+GENES_N$V10[l])/2, -40+Y_N, col="black") # V11 V10 start end crd
    points((GENES_N$V11[l]+GENES_N$V10[l])/2, -40+Y_N, col="blue", pch=20, cex=0.5)}}
  for (g in 1:nrow(GENES_N)) {
    X = c(ifelse(GENES_N$V5[g] == "+", GENES_N$start[g]-2e4, GENES_N$stop[g]+2e4), ifelse(GENES_N$V5[g] == "+", GENES_N$start[g]-2e4, GENES_N$stop[g]+2e4), ifelse(GENES_N$V5[g] == "+", GENES_N$start[g]+8e4, GENES_N$stop[g]-8e4))
    polygon(X, c(-90+Y_N, -70+Y_N, -80+Y_N), col="red", border=NA)}
  
  # CRD TCELL
  points(idx2pos1((CORR2_T$V1+CORR2_T$V2)/2), CORR2_T$V2-CORR2_T$V1+Y_T, col=rgb(0,0,1,abs(CORR2_T$V7)^2), pch=18, cex=0.6)
  for (m in 1:nrow(MODSs_T)) {
    polygon(c(idx2pos1(MODSs_T$LIDX[m]), idx2pos1((MODSs_T$LIDX[m]+MODSs_T$RIDX[m])/2), idx2pos1(MODSs_T$RIDX[m])), c(0+Y_T,MODSs_T$RIDX[m]-MODSs_T$LIDX[m]+Y_T, 0+Y_T), border="black", lwd=2)}
  for (m in 1:nrow(MODSs_T)) {
    polygon(c(idx2pos1(MODSs_T$LIDX[m]), idx2pos1(MODSs_T$RIDX[m]), MODSs_T$END[m], MODSs_T$START[m]), c(0+Y_T,0+Y_T,-40+Y_T,-40+Y_T), border=NA, col=rgb(0,0,1, alpha=0.3))
    arrows(MODSs_T$START[m], -40+Y_T, MODSs_T$END[m], -40+Y_T, length=0.025, angle=90, code=3, col="blue")}
  for (l in 1:nrow(GENES_T)) {if(GENES_T$associated[l]==1){
    segments(GENES_T$start[l], -80+Y_T, (GENES_T$V11[l]+GENES_T$V10[l])/2, -40+Y_T, col="black") # V11 V10 start end crd
    points((GENES_T$V11[l]+GENES_T$V10[l])/2, -40, col="blue", pch=20, cex=0.5)}}
  for (g in 1:nrow(GENES_T)) {
    X = c(ifelse(GENES_T$V5[g] == "+", GENES_T$start[g]-2e4, GENES_T$stop[g]+2e4), ifelse(GENES_T$V5[g] == "+", GENES_T$start[g]-2e4, GENES_T$stop[g]+2e4), ifelse(GENES_T$V5[g] == "+", GENES_T$start[g]+8e4, GENES_T$stop[g]-8e4))
    polygon(X, c(-90+Y_T, -70+Y_T, -80+Y_T), col="red", border=NA)}
  
  dev.off()
  
}

