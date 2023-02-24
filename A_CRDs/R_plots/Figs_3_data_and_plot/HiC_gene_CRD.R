#############################################################################################
#
# PLOT PCHiC support GENE CRD ASSOCIATIONS
#
#############################################################################################
# ~  /Users/dianaavalos/Programming/A_CRD_plots/HiCplots/HiC_gene_CRD.R

# Clean environment ------------------------------------
rm(list=ls())
gc()

# Packages ---------------------------------------------

library(qvalue)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(gridExtra)
library(cowplot)


Paper=c("#00AFBB", "#E7B800","#FC4E07")
color_palette= Paper
data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
conditions = list('mean')


#############################################################################################
#
# Functions
#
#############################################################################################


# Functions  ---------------------------------

get_corr_genes_formated <- function(corr_genes,genelist){
  
  corr_genes$distance = corr_genes$V4-corr_genes$V3
  corr_genes = corr_genes[,c(3,4,5,6,8,9)]
  colnames(corr_genes)[1:5] = c("pos1","pos2","gene1","gene2","pval")
  corr_genes$pval.adj = p.adjust(corr_genes$pval,method='fdr')
  corr_genes$sameCRD = -1
  setkeyv(corr_genes,c("gene1","gene2"))
  
  corr_genes = corr_genes[corr_genes$gene1 %in% genelist,]
  corr_genes = corr_genes[corr_genes$gene2 %in% genelist,]
  corr_genes$rowid = c(1:nrow(corr_genes))
  corr_genes
}


# HIC Functions  ---------------------------------
array_CRD_genes=mapdata_signif
compute_hic_validated <- function(PCHiC, array_CRD_genes){
  
  colnames(PCHiC)[1] = "baitChr"
  
  baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
  oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))
  
  genebed <- GRanges(seqnames=array_CRD_genes$phenotype_ID_chr,ranges=IRanges(start=array_CRD_genes$phenotype_ID_start, end=array_CRD_genes$phenotype_ID_end))
  CRDbed <- GRanges(seqnames=array_CRD_genes$CRD_ID_chr,ranges=IRanges(start=array_CRD_genes$CRD_ID_start, end=array_CRD_genes$CRD_ID_end))
  
  #fwd
  x = findOverlaps(baitbed,genebed)
  y = findOverlaps(oebed,CRDbed)
  tmp = rbind(as.data.frame(x),as.data.frame(y))
  validated.fwd = tmp[which(duplicated(tmp)),]
  
  #bwd
  x = findOverlaps(baitbed,CRDbed)
  y = findOverlaps(oebed,genebed)
  tmp = rbind(as.data.frame(x),as.data.frame(y)) 
  validated.bwd = tmp[which(duplicated(tmp)),]
  
  validated = unique(rbind(validated.fwd,validated.bwd))
  
  validated
}


compute_HiC_column_in_mapdata <- function(PCHiC, mapdata_signif,validated,cell_type){
  
  # select the HiC column of the cell type of interest
  currentHiCScore_all <- PCHiC[,c("Neu","Mon","nCD4")]
  names(currentHiCScore_all) <- c("neut", "mono", "tcell")
  currentHiCScore=currentHiCScore_all[[cell_type]]
  
  mapdata_signif_validated = rep(0,nrow(mapdata_signif))
  for(i in unique(validated$subjectHits)){
    currenthic = validated[validated$subjectHits==i,]$queryHits
    meanHiC=mean(currentHiCScore[currenthic])
    mapdata_signif_validated[i] = meanHiC
  }
  
  mapdata_signif_validated
}

 # G function
 # hic_validated = rep(1,nrow(PCHiC))
 # 
 # for(i in 1:nrow(validated)){
 #   currenthic = validated$queryHits[i]
 #   currentpval = mapdata_signif$bwd_pval[validated$subjectHits[i]]
 #   if(currentpval<hic_validated[currenthic]){
 #     hic_validated[currenthic] = currentpval
 #   }
 # }
 # 
 
 # #  D function but wrong
 #   mapdata_signif_validated = rep(0,nrow(mapdata_signif))
 #   for(i in 1:nrow(validated)){
 #     currenthic = validated$queryHits[i]
 #     currentmap = validated$subjectHits[i]
 #     
 #     currentHiCScore_all <- c(PCHiC$Neu[currenthic], PCHiC$Mon[currenthic], PCHiC$nCD4[currenthic])
 #     names(currentHiCScore_all) <- c("neut", "mono", "tcell")
 #     currentHiCScore=currentHiCScore_all[[cell_type]]
 #     
 #     if(currentHiCScore>mapdata_signif_validated[currentmap]){
 #       mapdata_signif_validated[currentmap] = currentHiCScore
 #     }
 #   }
 
 # D new
 # problem, old was overwriting the last HiC value for the CRD-gene overlap
 # better to take the average
 # loop over subjecthits values
 # take all the query hits associated
 # mean of query hits HiS


compute_ratio_hic_mapdata <-function(mapdata,CRDmindist,CRDmaxdist,cutoff=5){
  
  up = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist & mapdata$HIC>cutoff,na.rm=T)
  down = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist,na.rm=T)
  up/down*100
  
}

#############################################################################################
#
# DIRECTORIES AND FILES
#
#############################################################################################

# getting genelist
path_ref='/Users/dianaavalos/Programming/reference_files/'
protein_coding_genes = scan(paste0(path_ref, "gencode.v15.annotation.protein_coding.gene_id.txt"),what="")
long_nc_RNA_genes = scan(paste0(path_ref, "gencode.v15.annotation.long_noncoding_RNAs.gene_id.txt"),what="")
genelist = c(protein_coding_genes,long_nc_RNA_genes)

# paths 
path='/Users/dianaavalos/Programming/A_CRD_plots/RNA' 
path_out = '/Users/dianaavalos/Programming/A_CRD_plots/figs_geneCRD/'
path_CRD_gene_signifs='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/merged_TH/'
path_CRD_gene_notsignifs='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/nominal_merged/'

path_CRD='/Users/dianaavalos/Programming/A_CRD_plots/quantify_ALL/'
rna_file <- c('/EGAD00001002675_RNA.ALL.txt.gz', '/EGAD00001002674_RNA.ALL.txt.gz', '/EGAD00001002671_RNA.ALL.txt.gz')
names(rna_file) <- c("neut", "mono", "tcell")
PCHiC = fread('/Users/dianaavalos/Programming/HiC_nov20/PCHiC_peak_matrix_cutoff5.tsv')

path_out = '/Users/dianaavalos/Programming/A_CRD_plots/figs_geneCRD/'
path_out="/Users/dianaavalos/Programming/A_CRD_plots/HiCplots/"

data_type='hist'
cell_type='neut'
condition='mean'
FDR='0.05'

data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
condition = 'mean'

#############################################################################################
#
# MAIN LOOP
#
#############################################################################################

toplot_ALL=as.data.frame(matrix(NA, nrow = 0, ncol = 5))
Wilcoxsignif_ALL=as.data.frame(matrix(NA, nrow = 0, ncol = 3))


for(data_type in data_types){
  
  #TODO: create mapdatasignif for all cells and one datatype
  
  mapdata_signif_allcells=as.data.frame(matrix(NA, nrow = 0, ncol = 24))
  mapdata_notsignif_allcells=as.data.frame(matrix(NA, nrow = 0, ncol = 17))
  
  
  for(cell_type in cell_types){

      # names of files
      name_condition=paste0(data_type,'_',cell_type ,'_',condition)
      print(name_condition)
      file_CRD=paste0(path_CRD, data_type,'_',cell_type ,'.ALLchr.',condition,'.txt.gz')
      file_mapdata_signif=paste0(path_CRD_gene_signifs,data_type,'_',cell_type ,'_',condition,'_conditional.txt.gz')
      file_mapdata_notsignif=paste0(path_CRD_gene_notsignifs,data_type,'_',cell_type ,'_',condition,"_mapping_CRD_gene_ALL.txt.gz")
      
      # load files SIGNIFICANT GENE_CRD INTERACTION
      mapdata_signif = read.table(file_mapdata_signif, head=F, stringsAsFactors=F)
      
      colnames(mapdata_signif) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                   "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","rank",
                                   "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
      #Filter protein coding and long nc RNA genes
      mapdata_signif = mapdata_signif[mapdata_signif$phenotype_ID %in% genelist,]
      
      # load files NON SIGNIFICANT GENE_CRD INTERACTION
      mapdata_notsignif = read.table(file_mapdata_notsignif, head=F, stringsAsFactors=F)
      mapdata_notsignif = mapdata_notsignif[mapdata_notsignif$V1 %in% genelist,]
      colnames(mapdata_notsignif) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                   "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end",
                                   "bwd_pval","bwd_slope","flag")
      mapdata_notsignif=mapdata_notsignif[mapdata_notsignif$bwd_pval>=0.05,]
      
      
      ##### MAIN
      validated_signif=compute_hic_validated(PCHiC, mapdata_signif)
      validated_notsignif=compute_hic_validated(PCHiC, mapdata_notsignif)
      validated_signif=validated_signif
      validated_notsignif=validated_notsignif
      
      mapdata_signif$HIC=compute_HiC_column_in_mapdata(PCHiC, mapdata_signif,validated_signif,cell_type) # cell type for Hic
      mapdata_notsignif$HIC=compute_HiC_column_in_mapdata(PCHiC, mapdata_notsignif,validated_notsignif,cell_type)
      
      crd_dist_bins = c(0,1,1e04,2e04,5e04,1e05,2e05,5e05,1e06)
      mat_hic_signif = matrix(0,nrow=(length(crd_dist_bins)-1),ncol=1)
      mat_hic_notsignif = matrix(0,nrow=(length(crd_dist_bins)-1),ncol=1)
      
      for(i in 1:(length(crd_dist_bins)-1)){
        mat_hic_signif[i,1] = compute_ratio_hic_mapdata(mapdata_signif,crd_dist_bins[i],crd_dist_bins[i+1])
        mat_hic_notsignif[i,1] = compute_ratio_hic_mapdata(mapdata_notsignif,crd_dist_bins[i],crd_dist_bins[i+1])
      }
      
      
      mapdata_signif$cell=cell_type
      mapdata_notsignif$cell=cell_type
      mapdata_signif$data_type=data_type
      mapdata_notsignif$data_type=data_type
      
      mapdata_signif_allcells=rbind(mapdata_signif_allcells,mapdata_signif)
      mapdata_notsignif_allcells=rbind(mapdata_notsignif_allcells,mapdata_notsignif)
      
      #create summary table
      toplot_old = data.frame(hic_signif = mat_hic_signif[,1],hic_notsignif = mat_hic_notsignif[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
      toplot_1 = data.frame(hic = mat_hic_signif[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
      toplot_1$signif=1
      toplot_2 = data.frame(hic = mat_hic_notsignif[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
      toplot_2$signif=0
      toplot=rbind(toplot_1,toplot_2)
      toplot$cell=cell_type
      toplot$data_type=data_type
      
      toplot_ALL=rbind(toplot_ALL,toplot)
    }
  
   # for all data types now
    # SIGNIFICANCE OF THE DIFFERENCE TEST wilcox.test
    wilcox = data.frame(dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
    wilcox$pvalue=NA
    wilcox$data_type=data_type
    
      for(i in 1:(length(crd_dist_bins)-1)){
        print(i)
        mapdata_signif_distinterval = mapdata_signif_allcells[mapdata_signif_allcells$distance >= crd_dist_bins[i] & mapdata_signif_allcells$distance<crd_dist_bins[i+1],]
        mapdata_notsignif_distinterval = mapdata_notsignif_allcells[mapdata_notsignif_allcells$distance >= crd_dist_bins[i] & mapdata_notsignif_allcells$distance<crd_dist_bins[i+1],]
        TEST=wilcox.test(mapdata_signif_distinterval$HIC, mapdata_notsignif_distinterval$HIC)
        wilcox$pvalue[i]=TEST$p.value
      }
  
    Wilcoxsignif_ALL=rbind(Wilcoxsignif_ALL,wilcox)
    
}

Wilcoxsignif_ALL$signif=0
Wilcoxsignif_ALL[Wilcoxsignif_ALL$pvalue<0.05,]$signif=1  

# save values of toplot_ALL and Wilcoxsignif_ALL
write.table(Wilcoxsignif_ALL, file = paste0(path_out,"/Wilcoxsignif_ALL.txt"), sep = "\t",row.names = TRUE, col.names = NA, quote = F)
write.table(toplot_ALL, file = paste0(path_out,"/toplot_ALL.txt"), sep = "\t",row.names = TRUE, col.names = NA, quote = F)




# TO TAKE IT HALFWAY
toplot_ALL=read.table(file = paste0(path_out,"toplot_ALL.txt"), sep = "\t",header=T)
Wilcoxsignif_ALL=read.table(file = paste0(path_out,"Wilcoxsignif_ALL.txt"), sep = "\t",header=T)

# df3
# hCRD_s    mCRD_s   hCRD_ns   mCRD_ns      dist
# 1 35.480226 19.735615 42.285961 18.230295    inside
# 2 27.693659 15.747535 24.309987 12.824312    0-10kb
# 3 24.069487 12.792019 22.604362 11.217968   10-20kb
# 4 25.699019 16.252140 21.510054 11.529094   20-50kb
# 5 28.345250 18.480856 22.435096 13.123343  50-100kb
# 6 29.006622 16.665127 19.651058 13.332162 100-200kb
# 7 19.611260 11.268455 11.223430  7.288979 200-500kb
# 8  9.115591  7.610431  3.629758  2.220465   0.5-1Mb

df3=as.data.frame(matrix(NA, nrow = 8, ncol = 5))
colnames(df3)=c("hCRD_s","mCRD_s" , "hCRD_ns" ,"mCRD_ns", "dist")   
df3$dist= c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb")

hCRD_s=toplot_ALL %>% filter(data_type=="hist"& signif==1)
mCRD_s=toplot_ALL %>% filter(data_type=="methyl"& signif==1)
hCRD_ns=toplot_ALL %>% filter(data_type=="hist"& signif==0)
mCRD_ns=toplot_ALL %>% filter(data_type=="methyl"& signif==0)

for (d in c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb")){
  df3[df3$dist==d,]$hCRD_s =mean(hCRD_s[hCRD_s$dist==d,]$hic)
  df3[df3$dist==d,]$mCRD_s =mean(mCRD_s[mCRD_s$dist==d,]$hic)
  df3[df3$dist==d,]$hCRD_ns =mean(hCRD_ns[hCRD_ns$dist==d,]$hic)
  df3[df3$dist==d,]$mCRD_ns =mean(mCRD_ns[mCRD_ns$dist==d,]$hic)
}

#############################################################################################
#
# PLOT
#
#############################################################################################



df3 <- data.frame(hCRD_s=toplotH_S$mean_hCRD,mCRD_s=toplotM_S$mean_mCRD,hCRD_ns=toplotH_NS$mean_hCRD,mCRD_ns=toplotM_NS$mean_mCRD, dist=toplotM_S$dist)
df3$dist = factor(df3$dist,levels = df3$dist)
df.molten = reshape2::melt(df3,id.vars="dist")
df.molten$signif=1
df.molten$signif[which(df.molten$variable=="hCRD_s")]=0
df.molten$signif[which(df.molten$variable=="mCRD_s")]=0
df.molten$CRD_type="hCRD"
df.molten$CRD_type[which(df.molten$variable=="mCRD_s")]="mCRD"
df.molten$CRD_type[which(df.molten$variable=="mCRD_ns")]="mCRD"


df.molten$dist = factor(df.molten$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))

#############################################################################################
#
# PLOT with all columns
#
#############################################################################################



# other plot LAST NEW ONE UP AND DOWN
hCRDmolten <- df.molten %>% filter(CRD_type=="hCRD")
mCRDmolten <- df.molten %>% filter(CRD_type=="mCRD")
hCRDmolten$pvalue=rep(Wilcoxsignif_ALL[Wilcoxsignif_ALL$data_type=="hist",]$pvalue,2)
mCRDmolten$pvalue=rep(Wilcoxsignif_ALL[Wilcoxsignif_ALL$data_type=="methyl",]$pvalue,2)
mCRDmoltenS <- mCRDmolten %>% filter(variable=="mCRD_s")
hCRDmoltenS <- hCRDmolten %>% filter(variable=="hCRD_s")
for ( i in seq(1,length(df3$hCRD_s))){
  hCRDmoltenS$value[i]=max(df3$hCRD_s[i],df3$hCRD_ns[i])
}


a <- ggplot(hCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ ggtitle("PCHiC support for gene-CRD associations") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#00AFBB","#00AFBB"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ ylim(0,25) +
  geom_text(data = hCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+1,
              label = format(pvalue, scientific = T,digits = 2)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

b <- ggplot(mCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ 
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#E7B800","#E7B800"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ scale_y_reverse(limits=c(25,0))+
  geom_text(data = mCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+0.5,
              label = format(pvalue, scientific = T,digits = 2)))

# changing ylim to have big mCRD
c <- 
  ggplot(mCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ 
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#E7B800","#E7B800"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ scale_y_reverse(limits=c(15,0)) +
  geom_text(data = mCRDmoltenS,
                mapping = aes(
                  x = dist,
                  y = value+0.5,
                  label = format(pvalue, scientific = T,digits = 2)))


plot_grid(a, c, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/CRD_gene_HIC",".tiff"), units="in", width=10, height=10, dpi=300, compression = 'lzw')
plot_grid(a, b, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/CRD_gene_HIC_otherscale",".tiff"), units="in", width=10, height=10, dpi=300, compression = 'lzw')


#############################################################################################
#
# PLOT with some columns
#
#############################################################################################

df.molten <- df.molten %>%  filter(dist %in% c("50-100kb","100-200kb","200-500kb","0.5-1Mb"))
df3 <- df3 %>%  filter(dist %in% c("50-100kb","100-200kb","200-500kb","0.5-1Mb"))
Wilcoxsignif_ALL2 <- Wilcoxsignif_ALL %>%  filter(dist %in% c("50-100kb","100-200kb","200-500kb","0.5-1Mb"))

# other plot LAST NEW ONE UP AND DOWN
hCRDmolten <- df.molten %>% filter(CRD_type=="hCRD")
mCRDmolten <- df.molten %>% filter(CRD_type=="mCRD")
hCRDmolten$pvalue=rep(Wilcoxsignif_ALL2[Wilcoxsignif_ALL2$data_type=="hist",]$pvalue,2)
mCRDmolten$pvalue=rep(Wilcoxsignif_ALL2[Wilcoxsignif_ALL2$data_type=="methyl",]$pvalue,2)
mCRDmoltenS <- mCRDmolten %>% filter(variable=="mCRD_s")
hCRDmoltenS <- hCRDmolten %>% filter(variable=="hCRD_s")

for ( i in seq(1,length(df3$hCRD_s))){
  hCRDmoltenS$value[i]=max(df3$hCRD_s[i],df3$hCRD_ns[i])
}


a <- ggplot(hCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ ggtitle("PCHiC support for gene-CRD associations") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#00AFBB","#00AFBB"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ ylim(0,25) +
  geom_text(data = hCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+1,
              label = format(pvalue, scientific = T,digits = 2)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

b <- ggplot(mCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ 
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#E7B800","#E7B800"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ scale_y_reverse(limits=c(25,0))+
  geom_text(data = mCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+0.5,
              label = format(pvalue, scientific = T,digits = 2)))

# changing ylim
c <- ggplot(mCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ 
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#E7B800","#E7B800"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ scale_y_reverse(limits=c(15,0)) +
  geom_text(data = mCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+0.5,
              label = format(pvalue, scientific = T,digits = 2)))


plot_grid(a, c, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/CRD_gene_HIC","_fewcols.tiff"), units="in", width=10, height=10, dpi=300, compression = 'lzw')
plot_grid(a, b, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/CRD_gene_HIC_otherscale","_fewcols.tiff"), units="in", width=10, height=10, dpi=300, compression = 'lzw')



#############################################################################################
#
# PLOT cell types
#
#############################################################################################
hCRD_s=toplot_ALL %>% filter(data_type=="hist"& signif==1)
hCRD_s$dist = factor(hCRD_s$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))

ggplot(hCRD_s, aes(x = dist, y = hic,fill=factor(cell)))+ ggtitle("PCHiC support for gene-CRD associations") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =Paper)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ ylim(0,max(hCRD_s$hic)) 
  
  geom_text(data = hCRDmoltenS,
            mapping = aes(
              x = dist,
              y = hic+1,
              label = format(pvalue, scientific = T,digits = 2)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
