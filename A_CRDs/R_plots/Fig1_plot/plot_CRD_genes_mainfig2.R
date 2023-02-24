#############################################################################################
#
# PLOT gene CRD contacts
#
#############################################################################################

# Clean environment ------------------------------------

rm(list=ls())
gc()

# Packages ---------------------------------------------
library(biomaRt)
library(qvalue)
library(RColorBrewer)
library(data.table)
library(R.utils) 
library(dplyr)
library(ggplot2)
library(gggenes)
library("ggrepel")

# peak1000_in_bp=39423603
# peak2000_in_bp=60731675
# peak980_in_bp=39145760
# peak2020_in_bp=61606270

# Paths ---------------------------------------------

wdir="/Users/dianaavalos/Dropbox/CRD_immune/plot_CRD_genes_main_fig2"
setwd(wdir)
path=paste0(wdir,'/CRD')
path_peak=paste0(wdir,'/CORR')
path_associations=paste0(wdir,"/CRD_GENE")
outpath=paste0(wdir,"/OUT")
path_RNA=paste0(wdir,"/RNA")
rna_file <- c('EGAD00001002675_chr5.txt', 'EGAD00001002674_chr5.txt', 'EGAD00001002671_chr5.txt')
names(rna_file) <- c("neut", "mono", "tcell")


#############################################################################################
#
# Functions
#
#############################################################################################


# get gene info, name, position, strand, etc...
get_genes_info <- function(human, genelist) {
  gene_coords = biomaRt::getBM(
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id",
      "start_position",
      "end_position",'strand'
    ),
    filters = "ensembl_gene_id",
    values = genelist,
    mart = mart,
    useCache = FALSE
  )
  gene_coords$size = gene_coords$end_position - gene_coords$start_position
  gene_coords
}

# split genes to have the ENSEMBL name
get_genes_ENSG_ID <- function(array) {
  c <- strsplit(array, '\\.')
  Genes <- unlist(c)[2*(1:length(c))-1]
  Genes
}

# peaks to  bp for each peak of the CRD
mapping_peaks_to_bp_fct <- function (MODS, start=800, stop=2200) {
  array_peak_bp=as.data.frame(matrix(NA, nrow = 0, ncol = 2))
  names(array_peak_bp)<-c("peak","bp")
  for (peak in seq(start,stop)){
    array_peak_bp[nrow(array_peak_bp) + 1,] <- c(peak,unique(MODS[MODS$LIDX==peak,]$START))
  }
  array_peak_bp
}


#  extract for chrRNAseq info and mean gene activity
chr5_cell_genes_info <- function (path_RNA,rna_file,cell) {
  genes_chr5=fread(file.path(path_RNA,rna_file[[cell]]))
  genes_exp=rowSums(genes_chr5[,7:dim(genes_chr5)[2]])
  genes_chr5$meanexp=genes_exp
  genes_chr5$meanexplog=log10(genes_exp)
  genes_chr5_summary=cbind(genes_chr5[,c(1,2,3,4,5,6)],genes_chr5$meanexp,genes_chr5$meanexplog)
  colnames(genes_chr5_summary)=c("chr","start","stop","gene","info","strand","mean","logmean")
  genes_chr5_summary
}


# with TSS, find the closest peak
gene_TSS_to_closest_peak<- function (mapping_peaks_to_bp, genes_TSS=genes_chr5_inrange$start) {
  array_TSS_peak_bp=as.data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(array_TSS_peak_bp)<-c("TSS","min_dist","peak","bp")
  genes_peak=rep(-1,length(genes_TSS))
  for (idx_gene in seq(1,length(genes_TSS))){
    genes_peak[idx_gene]=mapping_peaks_to_bp[which.min(abs(mapping_peaks_to_bp$bp - genes_TSS[idx_gene])),]$peak
  }
  genes_peak
}
  
# Debug ---------------------------------------------

type='hist'
cell='neut'
chr_num=5
section=1000

#############################################################################################
#
# MAIN
#
#############################################################################################
# TODO: redo log10
# TODO: all in one plot
# TODO: genes associated, expressed and nothing

human <- useEnsembl(biomart = "ensembl", 
                            dataset = "hsapiens_gene_ensembl", 
                            mirror = "useast")

for (cell in c('neut','mono','tcell')){

  # for (chr_num in 1:22){#22
  
  # LOAD  ---------------------------------------------
  file = paste0(type, "_", cell, ".chr", chr_num, ".module.txt.gz")
  MODS = fread(file.path(path, file),
               head = TRUE,
               stringsAsFactors = FALSE)
  
  file_corr = paste0(type, "_", cell, ".corr.chr", chr_num, ".txt")
  CORR2 = fread(file.path(path_peak, file_corr),
                head = FALSE,
                stringsAsFactors = FALSE)
  
  file_assoc = paste0("FDR_0.05_",type,"_", cell,"_mean_mapping_CRD_gene_ALL.significant.txt")
  ASSOC = fread(file.path(path_associations, file_assoc),
                head = FALSE,
                stringsAsFactors = FALSE)
  
  
  # Data subsets  ---------------------------------------------
  
  
  # limit correlation for CIS and 250 sliding window
  CORR2$DIFF=CORR2$V2-CORR2$V1
  CORR2 <- subset(CORR2, DIFF <= 250)
  max_x=as.integer(max(CORR2$V2)/1000)*1000
  
  # I have a loop to split the plots in section of 1000 peaks disabled now
  # for (section in seq(1,max_x-1000,1000)) {
  print(paste0(type,"_",cell,"_chr",chr_num,' section ',section, ' ', section+1000))
  
  # MODS subsets  ---------------------------------------------
  
  CRD_count=5 # this is a filter on the CRD nodes
  # for (CRD_count in c(4,5)){
    MODS1 = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > CRD_count, ] # otherwise too many triangles
    MODS2 <- MODS1 %>% filter(RIDX > section & LIDX < section+1000) # create a subset of MODS in the window 1000-2000
    MODSs=MODS2
    CRDs_in_plot <- unique(MODSs$UUID)
    
    #  extract chr 5 RNAseq info and mean gene activity   ---------------------------------------------
    genes_chr5_summary=chr5_cell_genes_info(path_RNA,rna_file,cell)
    range_peaks=c( start=980, stop=2020 )
    mapping_peaks_to_bp=mapping_peaks_to_bp_fct(MODS, start=800, stop=2200)
    range_genes=mapping_peaks_to_bp %>% filter(peak %in% range_peaks)
    range_genes=range_genes$bp
    
    # SELECT GENES EXPRESSED AND ADD INFO BIOMART  ---------------------------------------------
    
    # select genes in the interval studied ( expressed)
    idx=genes_chr5_summary$start >= mapping_peaks_to_bp[1,]$bp &  genes_chr5_summary$start <= mapping_peaks_to_bp[dim(mapping_peaks_to_bp)[1],]$bp
    genes_chr5_inrange <- genes_chr5_summary[idx,] # 128 for neut
    
    # change gene TSS with closest peak and add mean activity
    genes_chr5_inrange$peak_TSS=gene_TSS_to_closest_peak(mapping_peaks_to_bp, genes_TSS=genes_chr5_inrange$start)
    # add info about genes
    genes_chr5_inrange$ensembl_gene_id <- get_genes_ENSG_ID( genes_chr5_inrange$gene)
    gene_coords = get_genes_info(human, genes_chr5_inrange$ensembl_gene_id)
    colnames(gene_coords)
    
    genes_chr5_withinfo=left_join(genes_chr5_inrange,gene_coords,by="ensembl_gene_id")
    
    # FILTER GENES
    
    # remove genes without hgnc_symbol
    genes_chr5_filtered <-  genes_chr5_withinfo[!(is.na(genes_chr5_withinfo$hgnc_symbol) | genes_chr5_withinfo$hgnc_symbol==""), ]
    summary(genes_chr5_filtered) # 99 genes in neut are expressed in our range
    
    ASSOC2 <- ASSOC %>% filter(V8 %in% CRDs_in_plot)  # 13 genes in neut have CRD-associations with the CRDs in the plot ( WE ONLY PLOT WITH COUNT=5)
    ASSOC %>% filter (V3 >= range_genes[1] & V3 <= range_genes[2]) # 1105 gene CRD associations in our range
    genes_chr5_filtered_associated <- inner_join(genes_chr5_filtered,ASSOC, by=c("gene"="V1")) # 27 genes are expressed and associated to a CRD (5%)
    genes_chr5_filtered_associated_withCRD_in_plot <- inner_join(genes_chr5_filtered,ASSOC2, by=c("gene"="V1")) # 7 genes are expressed and have CRD-associations with the CRDs in the plot 
    
    # summary for neut : 99 genes expressed and 27 associated and expressed and 13 associated with the CRDs in the plot because count=5, and 7 associated with a CRD in the plot + expressed
    # ONLY THE IMPORTANT ONES
    genes_chr5_filtered_associated_withCRD <- inner_join(genes_chr5_filtered,ASSOC, by=c("gene"="V1"))
    # 27 genes are expressed and associated to a CRD (5%) 
    # 99 genes are expressed in our range
    
    # THE GENES USED FOR THE PLOT
    # to filter for minimal mean expression?
    expressed_genes_100filter=genes_chr5_filtered_associated %>% filter(mean>=0)
    associated_genes_100filter=genes_chr5_filtered_associated_withCRD_in_plot %>% filter(mean>=0)
  

    #############################################################################################
    #
    # PLOT genes
    #
    #############################################################################################
    
    idx=!expressed_genes_100filter$gene %in% associated_genes_100filter$gene
    ALL_genes_100filter=expressed_genes_100filter[idx,]
    ALL_genes_100filter$associated=0
    ALL_genes_100filter$color_palette="azure4"
    associated_genes_100filter$associated=1
    associated_genes_100filter$color_palette="#cfa500"
    ALL_genes_100filter=rbind(ALL_genes_100filter,associated_genes_100filter)
    
    
    # array to use for plor
    ALL_genes_100filter$peak_TSS2=ALL_genes_100filter$peak_TSS+delta
    ALL_genes_100filter$zero=0
    # TODO: plot expressed genes and not only significantly associated genes!!!!
    delta=5 # width rectangles
    
    write.table(ALL_genes_100filter, file = paste0(outpath,"/txt/ALL_genes_100filter_","_",CRD_count,type,"_",cell,".chr",chr_num,"_",section,"2.txt"), sep = "\t",
                row.names = TRUE, col.names = NA)
    
    
}

CRD_count=5
for (cell in c('neut','mono','tcell')){
  ALL_genes_100filter=read.table(paste0("/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots","/txt/ALL_genes_100filter_","_",5,"hist","_",cell,".chr",chr_num,"_",section,"2.txt"), sep = "\t", header=T)
  
    
  ggplot() + 
    geom_rect(
      data = ALL_genes_100filter,
      mapping = aes(
        xmin = peak_TSS,
        xmax = peak_TSS2,
        ymin = zero,
        ymax = -logmean,
        fill = factor(associated)
      ),
      color = "azure2",
      alpha = 1
    ) +
    theme_classic() + scale_fill_manual(values = c("grey", "#E7B800")) +
      ylim(-12,0) +
      geom_text_repel(data = ALL_genes_100filter,
                      mapping = aes(
                        x = peak_TSS,
                        y = -logmean,
                        label = hgnc_symbol, angle=90
                      ), force_pull   = 0, # do not pull toward data points
                      nudge_y      = -1,
                      angle        = 90,
                      hjust        = 0,
                      max.overlaps = Inf,
                      color=ALL_genes_100filter$color_palette,
                      segment.size = 0.2,
                      max.iter = 1e4, max.time = 1
      ) +
      theme(
        axis.line.y  = element_blank(),
        axis.line.x  = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.text.y  = element_blank(),
        axis.title.y = element_blank()
      ) + scale_x_continuous(name = "x", breaks = seq(1000, 2000, 200)) +  scale_y_continuous(name ="y", breaks = -seq(1, 11, 2), labels=10^seq(1, 11, 2))
  
    ggsave(paste0(outpath,"/",type,"_", cell,"_",CRD_count,"_genes_expressed_and_associated",".tiff"), units="in", width=20, height=4, dpi=300, compression = 'lzw')
    ggsave(paste0(outpath,"/",type,"_", cell,"_",CRD_count,"_genes_expressed_and_associated2",".tiff"), units="in", width=20, height=2, dpi=300, compression = 'lzw')
}

    #############################################################################################
    #
    # PLOT CRDS
    #
    #############################################################################################
    
    # add MODS
    
      #CRDS_to_remove=c("5_internal_12901","5_internal_11109") # because they are one inside the other
      
      CORR2sub <- CORR2 %>% filter(V1 >=section -250 & V2<=section+1000+250)
      range <- c(min(CORR2sub$V3),max(CORR2sub$V4))
      MODS1 = MODS[MODS$N_REG > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ] # otherwise too many triangles
      MODSs <- MODS1 %>% filter(RIDX > section & LIDX < section+1000) # create a subset of MODS in the window 1000-2000
      CRDs_in_plot <- unique(MODSs$UUID)
      
      # PLOT CRD
      filename=paste0(type,"_",cell,".chr",chr_num,"_",section,"_",CRD_count,"_figCRD_names.png")
      png(paste0(outpath,"/",filename), width=3000, height = 1000, res = 200) #  type="cairo",with limits
      plot((CORR2sub$V2+CORR2sub$V1)/2, CORR2sub$V2-CORR2sub$V1, pch=18, col=rgb(0,0,1,abs(CORR2sub$V7)^2), xlab="corr", xlim=c(section, section +1000), ylim=c(0, 250)) #, bty="n", ylab="", xlab="Chromatin peak index on chromosome 14", yaxt="n")
      abline(h=0)
      for (m in 1:nrow(MODSs)) {
        polygon(c(MODSs$LIDX[m], (MODSs$LIDX[m]+MODSs$RIDX[m])/2, MODSs$RIDX[m]), c(0,MODSs$RIDX[m]-MODSs$LIDX[m], 0), border="black", lwd=2)
        text(x=(MODSs$LIDX[m]+MODSs$RIDX[m])/2, y=MODSs$RIDX[m]-MODSs$LIDX[m]+10,labels=paste0(MODSs$UUID[m],"_",m),col = "red")}
      dev.off()
      
      # idx=!MODSs$UUID%in% CRDS_to_remove
      # MODSs2=MODSs[idx,]
      MODSs2=MODSs
      filename=paste0(type,"_",cell,".chr",chr_num,"_",section,"_",CRD_count,"_figCRD.png")
      png(paste0(outpath,"/",filename), width=3000, height = 1000, res = 200) # with limits
      plot((CORR2sub$V2+CORR2sub$V1)/2, CORR2sub$V2-CORR2sub$V1, pch=18, col=rgb(0,0,1,abs(CORR2sub$V7)^2), xlab="corr", xlim=c(section, section +1000), ylim=c(0, 250)) #, bty="n", ylab="", xlab="Chromatin peak index on chromosome 14", yaxt="n")
      abline(h=0)
      for (m in 1:nrow(MODSs2)) {
        polygon(c(MODSs2$LIDX[m], (MODSs2$LIDX[m]+MODSs2$RIDX[m])/2, MODSs2$RIDX[m]), c(0,MODSs2$RIDX[m]-MODSs2$LIDX[m], 0), border="black", lwd=2)}
      dev.off()
  # }
}



# old scripts, print arrows
# careful we take beggining and end from biomart, not TSS

# min size of gene
idx=ALL_genes_100filter$length < 160000
ALL_genes_100filter$end_position_fake=ALL_genes_100filter$V4
ALL_genes_100filter[idx,]$end_position_fake=ALL_genes_100filter[idx,]$start_position + 160000

ALL_genes_100filter=ALL_genes_100filter[order(associated),]
# also not associated genes in same space
index=ALL_genes_100filter$associated==0
color_palette= rep("grey",length(ALL_genes_100filter$start_position))
color_palette[index]="#E69F00"
ALL_genes_100filter$color_palette=color_palette
ALL_genes_100filter$cell_type=cell

ggplot(
  ALL_genes_100filter,
  aes(xmin = start_position, xmax = end_position_fake, y = cell_type, fill = hgnc_symbol, label = hgnc_symbol))  +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(1, "mm")) +
  xlim(range_genes[1],range_genes[2]) +
  geom_blank(data = ALL_genes_100filter) +
  scale_fill_manual(values =ALL_genes_100filter$color_palette) +
  geom_text_repel(x=ALL_genes_100filter$start_position,aes(label = hgnc_symbol),size = 3,max.overlaps = Inf,angle = 90, force_pull   = 0, # do not pull toward data points
                  nudge_y      = -0.2,
                  hjust        = 0,
                  color=ALL_genes_100filter$color_palette,
                  segment.size = 0.05)+
  theme_genes()
#annotate(geom="text", x=ALL_genes_100filter$start_position, y=1.2, label=ALL_genes_100filter$hgnc_symbol,color="red",size = 2)   +

ggsave(paste0(outpath,"/",type,"_", cell,"_genes_with_CRDs_not_assoc3",".tiff"), units="in", width=20, height=4, dpi=300, compression = 'lzw')



