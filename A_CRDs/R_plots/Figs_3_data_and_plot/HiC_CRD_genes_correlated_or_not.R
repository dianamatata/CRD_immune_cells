
#############################################################################################
#
# PLOT PCHiC support GENEs COEXPRESSED OR NOT AND ASSOCIATED WITH CRDs
#
#############################################################################################
# ~  /Volumes/Elements/Programming\ PhD/A_CRD_plots/HiCplots/...

# Clean environment ------------------------------------
rm(list=ls())
gc()


# Packages ---------------------------------------------

library(qvalue)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(tidyverse)
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



compute_ratio_hic_mapdata <-function(mapdata,CRDmindist,CRDmaxdist,cutoff=5){
  
  up = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist & mapdata$HIC>cutoff,na.rm=T)
  down = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist,na.rm=T)
  up/down*100
  
}

#HiC data

mindist_btw_genes=0
maxdist_btw_genes=1e09
i=1
mindist_CRD_gene=gene_dist_bins[i]
maxdist_CRD_gene=gene_dist_bins[i+1]
pval.cutoff=0.01
threshold=5


compute_ratio_hic_2022 <-function(coexpressed=T,mindist_btw_genes,maxdist_btw_genes,mindist_CRD_gene,maxdist_CRD_gene,df_associations2,pval.cutoff=0.01,threshold=5){
  if(coexpressed){
    # up= genes coexpressed  from same CRD
    up = sum(
      df_associations2$meanCRD_dist_2_genes >= mindist_CRD_gene &
        df_associations2$meanCRD_dist_2_genes < maxdist_CRD_gene &
        abs(df_associations2$dist_btw_genes) >= mindist_btw_genes &
        abs(df_associations2$dist_btw_genes) < maxdist_btw_genes &
        df_associations2$adj_pval < pval.cutoff &
        df_associations2$meanHiC_CRD_genes > threshold ,
      na.rm = T
    )
    
    # down= genes coexpressed from same CRD but not Hic significant
    down = sum(
      df_associations2$meanCRD_dist_2_genes >= mindist_CRD_gene &
        df_associations2$meanCRD_dist_2_genes < maxdist_CRD_gene &
        abs(df_associations2$dist_btw_genes) >= mindist_btw_genes &
        abs(df_associations2$dist_btw_genes) < maxdist_btw_genes &
        df_associations2$adj_pval < pval.cutoff,
      na.rm = T
    )
    
    } else {
      
      up = sum(
      # pvalue superior. not coexpressed    up = sum(
      df_associations2$meanCRD_dist_2_genes >= mindist_CRD_gene &
        df_associations2$meanCRD_dist_2_genes < maxdist_CRD_gene &
        abs(df_associations2$dist_btw_genes) >= mindist_btw_genes &
        abs(df_associations2$dist_btw_genes) < maxdist_btw_genes &
        df_associations2$adj_pval >= pval.cutoff &
        df_associations2$meanHiC_CRD_genes > threshold,
      na.rm = T
      )
  
    # down= genes coexpressed from same CRD but not Hic significant
    down = sum(
      df_associations2$meanCRD_dist_2_genes >= mindist_CRD_gene &
        df_associations2$meanCRD_dist_2_genes < maxdist_CRD_gene &
        abs(df_associations2$dist_btw_genes) >= mindist_btw_genes &
        abs(df_associations2$dist_btw_genes) < maxdist_btw_genes &
        df_associations2$adj_pval  >=  pval.cutoff,
      na.rm = T
    )
  
    }
   up/down*100

}


#HiC data, after annotate_map_data_with_corr_genes_OLD

compute_ratio_hic_OLD <-function(corr_genes,coexpressed=T,CRDmindist,CRDmaxdist,mindist,maxdist,pval.cutoff=0.01,threshold=5){
  if(coexpressed){
    up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$hic>threshold & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    down = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
  } else {
    up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$hic>threshold & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    down = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
  }
  up/down*100
}  


# Annotate corr_genes with CRD mapping + HIC ---------------------------------
# new HiC

annotate_map_data_with_corr_genes_2022 <- function(mapdata, corr_genes){
  
  #Annotate corr_genes with CRD mapping
  CRD_IDs = unique(mapdata$CRD_ID)
  df_associations=data.frame(matrix(ncol = 15, nrow = 0))
  corr_genes=data.table(corr_genes)
  colnames(df_associations) <- c("CRD_ID", "CRD_ID_start","CRD_ID_end","gene1", "gene2","gene1_start", "gene2_start","adj_pval","dist_btw_genes","meanCRD_dist_2_genes","CRD_gene1_dist","CRD_gene2_dist", "pval_CRDgene1","pval_CRDgene2","meanHiC_CRD_genes")
  # loop over CRDs to identify the genes correlated to this CRD and also if these 2 genes are in the correlated gene pairs
  
  for(i in 1:length(CRD_IDs)){
    if(i%%100==0){print(i)}
    idx = which(mapdata$CRD_ID == CRD_IDs[i])
    if(length(idx)>1){
      for(k in 1:(length(idx)-1)){
        for(l in (k+1):length(idx)){
          geneA = mapdata$phenotype_ID[idx[k]]
          geneB = mapdata$phenotype_ID[idx[l]]
          idx_association_2corr_genes_1_CRD =  c(unlist(corr_genes[gene1 == geneB & gene2 == geneA,9]),unlist(corr_genes[gene1 == geneA & gene2 == geneB,9]))
          if(length(idx_association_2corr_genes_1_CRD)>0){ 
            # if they are co-expressed 
            pval_corr_genes=corr_genes[idx_association_2corr_genes_1_CRD,]$pval.adj
          } else {
            # if they are not but regulated by the same CRD
            pval_corr_genes=1
          }
          
          meanCRD_dist_2_genes = mean(c(abs(mapdata$distance_CRD_gene[idx[k]]),abs(mapdata$distance_CRD_gene[idx[l]])) )# distance of the gene to the CRD CRD_IDs[i]
          dist_btw_genes = abs(mapdata$phenotype_ID_start[idx[l]]-mapdata$phenotype_ID_start[idx[k]])# distance between the genes
          a = c(
            CRD_IDs[i],
            mapdata$CRD_ID_start[idx[k]],
            mapdata$CRD_ID_end[idx[k]],
            geneA,
            geneB,
            mapdata$phenotype_ID_start[idx[k]],
            mapdata$phenotype_ID_start[idx[l]],
            pval_corr_genes,
            dist_btw_genes,
            meanCRD_dist_2_genes,
            abs(mapdata$distance_CRD_gene[idx[k]]),
            abs(mapdata$distance_CRD_gene[idx[l]]),
            mapdata$bwd_pval[idx[k]],
            mapdata$bwd_pval[idx[l]],
            meanHiC_CRD_genes = mean(mapdata$HIC[idx[k]],mapdata$HIC[idx[l]])
          )
          df_associations[nrow(df_associations) + 1,]=a
        }
      }
    }
  }
  
  df_associations$dist_btw_genes=as.integer(df_associations$dist_btw_genes)
  df_associations$meanCRD_dist_2_genes=as.integer(df_associations$meanCRD_dist_2_genes)
  df_associations$sameCRD=0
  df_associations$gene1_start=as.integer(df_associations$gene1_start)
  df_associations$gene2_start=as.numeric(df_associations$gene2_start)
  df_associations$adj_pval=as.numeric(df_associations$adj_pval)
  df_associations
}

# Annotate corr_genes with CRD mapping + HIC ---------------------------------
# in corrgenes we overwrite half of the genes since they are associated with 2 CRDs, pb...
annotate_map_data_with_corr_genes_OLD <- function(mapdata_signif, corr_genes){
  
  #Annotate corr_genes with CRD mapping
  
  CRD_IDs = unique(mapdata_signif$CRD_ID)
  IDXs = c()
  CRDdistance = c()
  CRDdistance_gene1=c()
  CRDdistance_gene2=c()
  hicsupport = c()
  corr_genes$CRDdistance_gene1=-1
  corr_genes$CRDdistance_gene2=-1
  
  for(i in 1:length(CRD_IDs)){
    # cat(i,"\n")
    idx = which(mapdata_signif$CRD_ID == CRD_IDs[i])
    if(length(idx)>1){
      for(k in 1:(length(idx)-1)){
        for(l in (k+1):length(idx)){
          geneA = mapdata_signif$phenotype_ID[idx[k]]
          geneB = mapdata_signif$phenotype_ID[idx[l]]
          myhit = c(unlist(corr_genes[gene1 == geneB & gene2 == geneA,9]),unlist(corr_genes[gene1 == geneA & gene2 == geneB,9]))
          if(length(myhit)>0){
            IDXs = c(IDXs,myhit)
            tmp = mean(abs(mapdata_signif$distance[idx[k]]),abs(mapdata_signif$distance[idx[l]])) # distance of the gene to the CRD CRD_IDs[i]
            tmp2 = mean(mapdata_signif$HIC[idx[k]],mapdata_signif$HIC[idx[l]])
            hicsupport = c(hicsupport,tmp2)
            
            CRDdistance = c(CRDdistance,tmp) # mean distance with the 2 genes
            CRDdistance_gene1 = c(CRDdistance_gene1,abs(mapdata_signif$distance[idx[k]])) # dist of gene 1
            CRDdistance_gene2 = c(CRDdistance_gene2,abs(mapdata_signif$distance[idx[l]])) # dist of gene 2
          }
        }
      }
    }
  }
  
  # > length(unique(IDXs))
  # [1] 9236
  # > length(IDXs)
  # [1] 14989
  # in corrgenes we overwrite half of the genes since they are associated with 2 CRDs
  
  tmpdf = data.table(idx=IDXs,CRDdistance=CRDdistance,hicsupport=hicsupport, CRDdistance_gene1=CRDdistance_gene1, CRDdistance_gene2=CRDdistance_gene2)
  
  corr_genes$sameCRD[tmpdf$idx] = tmpdf$CRDdistance # by default corr_genes$sameCRD=-1. replace it with the min distance to the CRDs
  corr_genes$hic = 0
  corr_genes$hic[tmpdf$idx] = tmpdf$hicsupport
  corr_genes$CRDdistance_gene1[tmpdf$idx] = tmpdf$CRDdistance_gene1
  corr_genes$CRDdistance_gene2[tmpdf$idx] = tmpdf$CRDdistance_gene2
  corr_genes
}



#############################################################################################
#
# DIRECTORIES AND FILES
#
#############################################################################################

# getting genelist
path_ref='/Volumes/Elements/reference_files/'
protein_coding_genes = scan(paste0(path_ref, "gencode.v15.annotation.protein_coding.gene_id.txt"),what="")
long_nc_RNA_genes = scan(paste0(path_ref, "gencode.v15.annotation.long_noncoding_RNAs.gene_id.txt"),what="")
genelist = c(protein_coding_genes,long_nc_RNA_genes)

# paths 
path='/Volumes/Elements/Programming\ PhD/A_CRD_plots/RNA' 
path_out = '/Volumes/Elements/Programming\ PhD/A_CRD_plots/figs_geneCRD/'
path_CRD_gene_signifs='/Volumes/Elements/Programming\ PhD/A_CRD_plots/CRD_genes_5_for_plots/merged_TH/'
path_CRD_gene_notsignifs='/Volumes/Elements/Programming\ PhD/A_CRD_plots/CRD_genes_5_for_plots/nominal_merged/'

path_CRD='/Volumes/Elements/Programming\ PhD/A_CRD_plots/quantify_ALL/'
rna_file <- c('/EGAD00001002675_RNA.ALL.txt.gz', '/EGAD00001002674_RNA.ALL.txt.gz', '/EGAD00001002671_RNA.ALL.txt.gz')
names(rna_file) <- c("neut", "mono", "tcell")
PCHiC = fread('/Volumes/Elements/Programming\ PhD/HiC_nov20/PCHiC_peak_matrix_cutoff5.tsv')

path_out = '/Volumes/Elements/Programming\ PhD/A_CRD_plots/figs_geneCRD/'
path_out="/Volumes/Elements/Programming\ PhD/A_CRD_plots/HiCplots/"

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
  
  corr_genes2_allcells=as.data.frame(matrix(NA, nrow = 0, ncol = 24))
  
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
                                   "nb_variants","distance_CRD_gene","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","rank",
                                   "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
      
      # load corr genes
      corr_genes=fread(paste0(path,rna_file[[cell_type]]))
      corr_genes=get_corr_genes_formated(corr_genes,genelist)
      
      #Filter protein coding and long nc RNA genes
      mapdata_signif = mapdata_signif[mapdata_signif$phenotype_ID %in% genelist,]
      
  
      ##### MAIN
      validated_signif=compute_hic_validated(PCHiC, mapdata_signif)
      mapdata_signif$HIC=compute_HiC_column_in_mapdata(PCHiC, mapdata_signif,validated_signif,cell_type) # cell type for Hic
  
      # NOW we have it written for all cells
      df_associations2=annotate_map_data_with_corr_genes_2022(mapdata_signif, corr_genes)
      #write.table(df_associations,paste0("/Volumes/Elements/Programming\ PhD/A_CRD_plots/RNA/df_associations_",name_condition,".txt"),sep="\t",row.names=FALSE, quote=FALSE)
      # df_associations2=read.table(paste0("/Volumes/Elements/Programming\ PhD/A_CRD_plots/RNA/df_associations_",name_condition,".txt"),sep="\t", header=T)
      head(df_associations2)
      
      corr_genes2=annotate_map_data_with_corr_genes_OLD(mapdata_signif, corr_genes)
      # problem of annotate_map_data_with_corr_genes_OLD: we overwrite gene associations with 2 CRDs... but so far this is what we use..
      
  
      gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06,3e09)
      coexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)
      notcoexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)
      
      for(i in 1:(length(gene_dist_bins)-1)){
        coexpressed_mat_hic[i,1] = compute_ratio_hic_OLD(corr_genes2,T,0,1e09,gene_dist_bins[i],gene_dist_bins[i+1])
        notcoexpressed_mat_hic[i,1] =  compute_ratio_hic_OLD(corr_genes2,F,0,1e09,gene_dist_bins[i],gene_dist_bins[i+1])
      }
      
      
      # save all corr_genes to compute the mean later 
      corr_genes2$cell=cell_type
      corr_genes2$data_type=data_type
      corr_genes2_allcells=rbind(corr_genes2_allcells,corr_genes2)
      
      
      #create summary table
      toplot_old = data.frame(hic_signif = coexpressed_mat_hic[,1],hic_notsignif = notcoexpressed_mat_hic[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
      toplot_1 = data.frame(hic = coexpressed_mat_hic[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
      toplot_1$signif=1
      toplot_2 = data.frame(hic = notcoexpressed_mat_hic[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
      toplot_2$signif=0
      toplot=rbind(toplot_1,toplot_2)
      toplot$cell=cell_type
      toplot$data_type=data_type
      toplot_ALL=rbind(toplot_ALL,toplot)
    
   }
  
  # for all data types now
  # SIGNIFICANCE OF THE DIFFERENCE TEST wilcox.test
  wilcox = data.frame(dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
  wilcox$pvalue=NA
  wilcox$data_type=data_type
  pval.cutoff=0.01
  threshold=5
  
  gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06,3e09)
  
  for(i in 1:(length(gene_dist_bins)-1)){
    coexpressed_genes_distinterval=corr_genes2_allcells[corr_genes2_allcells$sameCRD!=-1 & # associated with sameCRD
                                                          corr_genes2_allcells$pval.adj < pval.cutoff  & # coexpressed
                                                          corr_genes2_allcells$distance >= gene_dist_bins[i] & # same distbin
                                                          corr_genes2_allcells$distance<gene_dist_bins[i+1],]
    
    notcoexpressed_genes_distinterval=corr_genes2_allcells[corr_genes2_allcells$sameCRD!=-1 & # associated with sameCRD
                                                          corr_genes2_allcells$pval.adj >= pval.cutoff  & # coexpressed
                                                          corr_genes2_allcells$distance >= gene_dist_bins[i] & # same distbin
                                                          corr_genes2_allcells$distance<gene_dist_bins[i+1],]
    TEST=wilcox.test(coexpressed_genes_distinterval$hic, notcoexpressed_genes_distinterval$hic)
    wilcox$pvalue[i]=TEST$p.value
  }
  
  Wilcoxsignif_ALL=rbind(Wilcoxsignif_ALL,wilcox)
  
}

Wilcoxsignif_ALL$signif=0
Wilcoxsignif_ALL[Wilcoxsignif_ALL$pvalue<0.05,]$signif=1  



# save values of toplot_ALL and Wilcoxsignif_ALL
write.table(Wilcoxsignif_ALL, file = paste0(path_out,"/Wilcoxsignif_ALL_coexpressedgenes.txt"), sep = "\t",row.names = TRUE, col.names = NA, quote = F)
write.table(toplot_ALL, file = paste0(path_out,"/toplot_ALL_coexpressedgenes.txt"), sep = "\t",row.names = TRUE, col.names = NA, quote = F)

Wilcoxsignif_ALL=read.table(file = paste0(path_out,"/Wilcoxsignif_ALL_coexpressedgenes.txt"), sep = "\t",header=T)
toplot_ALL=read.table(file = paste0(path_out,"/toplot_ALL_coexpressedgenes.txt"), sep = "\t",header=T)



df3=as.data.frame(matrix(NA, nrow = 8, ncol = 5))
colnames(df3)=c("hCRD_s","mCRD_s" , "hCRD_ns" ,"mCRD_ns", "dist")   
df3$dist= c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb")

hCRD_s=toplot_ALL %>% filter(data_type=="hist"& signif==1)
mCRD_s=toplot_ALL %>% filter(data_type=="methyl"& signif==1)
hCRD_ns=toplot_ALL %>% filter(data_type=="hist"& signif==0)
mCRD_ns=toplot_ALL %>% filter(data_type=="methyl"& signif==0)

for (d in c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb")){
  df3[df3$dist==d,]$hCRD_s =mean(hCRD_s[hCRD_s$dist==d,]$hic)
  df3[df3$dist==d,]$mCRD_s =mean(mCRD_s[mCRD_s$dist==d,]$hic)
  df3[df3$dist==d,]$hCRD_ns =mean(hCRD_ns[hCRD_ns$dist==d,]$hic)
  df3[df3$dist==d,]$mCRD_ns =mean(mCRD_ns[mCRD_ns$dist==d,]$hic)
}

# careful the dist is not updated yet
df3$dist=c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb")


df3$dist = factor(df3$dist,levels = df3$dist)

df.molten = reshape2::melt(df3,id.vars="dist")
df.molten$signif=1
df.molten$signif[which(df.molten$variable=="hCRD_s")]=0
df.molten$signif[which(df.molten$variable=="mCRD_s")]=0
df.molten$CRD_type="hCRD"
df.molten$CRD_type[which(df.molten$variable=="mCRD_s")]="mCRD"
df.molten$CRD_type[which(df.molten$variable=="mCRD_ns")]="mCRD"


df.molten$dist = factor(df.molten$dist,levels = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))


print("PCHiC support for gene pairs associated with the same CRD")

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
  mCRDmoltenS$value[i]=max(df3$mCRD_s[i],df3$mCRD_ns[i])
  
}


a <- ggplot(hCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ ggtitle("PCHiC support for gene-gene associations") +
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
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ scale_y_reverse(limits=c(18,0)) +
  geom_text(data = mCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+0.5,
              label = format(pvalue, scientific = T,digits = 2)))


plot_grid(a, c, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/gene_gene_HIC",".tiff"), units="in", width=10, height=10, dpi=300, compression = 'lzw')
plot_grid(a, b, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/gene_gene_HIC_otherscale",".tiff"), units="in", width=10, height=10, dpi=300, compression = 'lzw')


#############################################################################################
#
# PLOT with SOME columns
#
#############################################################################################
dist_subset=c("0.2-0.5Mb","0.5-1Mb",">1Mb")
df.molten <- df.molten %>%  filter(dist %in% dist_subset)
df3 <- df3 %>%  filter(dist %in% dist_subset)
Wilcoxsignif_ALL2 <- Wilcoxsignif_ALL %>%  filter(dist %in% dist_subset)

# other plot LAST NEW ONE UP AND DOWN
hCRDmolten <- df.molten %>% filter(CRD_type=="hCRD")
mCRDmolten <- df.molten %>% filter(CRD_type=="mCRD")
hCRDmolten$pvalue=rep(Wilcoxsignif_ALL2[Wilcoxsignif_ALL2$data_type=="hist",]$pvalue,2)
mCRDmolten$pvalue=rep(Wilcoxsignif_ALL2[Wilcoxsignif_ALL2$data_type=="methyl",]$pvalue,2)
mCRDmoltenS <- mCRDmolten %>% filter(variable=="mCRD_s")
hCRDmoltenS <- hCRDmolten %>% filter(variable=="hCRD_s")

for ( i in seq(1,length(df3$hCRD_s))){
  hCRDmoltenS$value[i]=max(df3$hCRD_s[i],df3$hCRD_ns[i])
  mCRDmoltenS$value[i]=max(df3$mCRD_s[i],df3$mCRD_ns[i])
  
}


a <- ggplot(hCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ ggtitle("PCHiC support for gene-gene associations") +
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
  labs(x = "CRD-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ scale_y_reverse(limits=c(18,0)) +
  geom_text(data = mCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+0.5,
              label = format(pvalue, scientific = T,digits = 2)))


plot_grid(a, c, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/gene_gene_HIC","_fewcols.tiff"), units="in", width=6, height=10, dpi=300, compression = 'lzw')
plot_grid(a, b, ncol = 1, nrow = 2)
ggsave(paste0(path_out,"/gene_gene_HIC_otherscale","_fewcols.tiff"), units="in", width=6, height=10, dpi=300, compression = 'lzw')











