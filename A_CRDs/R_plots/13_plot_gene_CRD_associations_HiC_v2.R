#**************************** PLOTS WITH  Hi-C analysis

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

Paper=c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07","#972D15")
color_palette= Paper


#############################################################################################
#
# FUNCTIONS
#
#############################################################################################

compute_ratio <-function(coexpressed=T,mindist,maxdist,CRDmindist,CRDmaxdist,corr_genes,pval.cutoff=0.01){
  if(coexpressed){
    up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    down = sum(corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
  } else {
    up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    down = sum(corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
  }
  up/down*100
}
# new 
compute_ratio_n <-function(coexpressed=T,mindist,maxdist,CRDmindist,CRDmaxdist,corr_genes,pval.cutoff=0.01){
  if(coexpressed){
    # corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist
    up_both = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 ==0) & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    up_one = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 !=0) & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    down = sum(corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
  } else {
    up_both = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 ==0) & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    up_one = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 !=0) & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
    down = sum(corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
  }
  c(up_both/down*100, up_one/down*100)
}

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


#############################################################################################
#
# INTEGRATION OF HIC WITH TRANS DATA by with compute_hic_validated
#
#############################################################################################

# slightly different from the one in 26_plot_trans_analysis.R

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


compute_HiC_column_in_mapdata <- function(PCHiC, mapdata,validated,cell_type){
  
  hic_validated = rep(1,nrow(PCHiC))
  
  for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentpval = mapdata$bwd_pval[validated$subjectHits[i]]
    if(currentpval<hic_validated[currenthic]){
      hic_validated[currenthic] = currentpval
    }
  }
  

  mapdata_validated = rep(0,nrow(mapdata))
  for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentmap = validated$subjectHits[i]
    
    currentHiCScore_all <- c(PCHiC$Neu[currenthic], PCHiC$Mon[currenthic], PCHiC$nCD4[currenthic])
    names(currentHiCScore_all) <- c("neut", "mono", "tcell")
    currentHiCScore=currentHiCScore_all[[cell_type]]
    
    if(currentHiCScore>mapdata_validated[currentmap]){
      mapdata_validated[currentmap] = currentHiCScore
    }
  }
  
  mapdata_validated
}


compute_ratio_hic_mapdata <-function(CRDmindist,CRDmaxdist,cutoff=5){
  
  up = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist & mapdata_validated>cutoff,na.rm=T)
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

path='/Users/dianaavalos/Programming/Hi-C_correlated_peaks/'
path_out = '/Users/dianaavalos/Programming/A_CRD_plots/figs_geneCRD/'
path_CRD_genes='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/merged_TH/'


path_CRD='/Users/dianaavalos/Programming/A_CRD_plots/quantify_ALL/'
rna_file <- c('/EGAD00001002675_RNA.ALL.txt.gz', '/EGAD00001002674_RNA.ALL.txt.gz', '/EGAD00001002671_RNA.ALL.txt.gz')
names(rna_file) <- c("neut", "mono", "tcell")

PCHiC = fread('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/PCHiC_peak_matrix_cutoff5.tsv')

data_type='hist'
cell_type='neut'
condition='mean'
FDR='0.05'


data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
conditions = list('mean', 'loom')
for(data_type in data_types){
  for(cell_type in cell_types){
    for(condition in conditions){
      
      # names of files
      name_condition=paste0(data_type,'_',cell_type ,'_',condition)
      file_CRD=paste0(path_CRD, data_type,'_',cell_type ,'.ALLchr.',condition,'.txt.gz')
      file_mapdata=paste0(path_CRD_genes,data_type,'_',cell_type ,'_',condition,'_conditional.txt.gz')
      print(name_condition)
      
      # download files
      allCRDs = fread(file_CRD,header=F)
      mapdata = read.table(file_mapdata, hea=F, stringsAsFactors=F)
      corr_genes=fread(paste0(path,rna_file[[cell_type]]))
      corr_genes=get_corr_genes_formated(corr_genes,genelist)
      
      colnames(mapdata) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                    "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","rank",
                                    "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
      
      
      #Filter protein coding and long nc RNA genes
      mapdata = mapdata[mapdata$phenotype_ID %in% genelist,]
      mapdata = mapdata[order(mapdata[,2],mapdata[,3]),]
      nb_CRD_not_associated = nrow(allCRDs) - length(unique(mapdata$CRD_ID))

      ##### MAIN
      nb_genes_not_associated = length(unique(corr_genes$gene1)) - length(unique(mapdata$phenotype_ID))
      validated=compute_hic_validated(PCHiC, mapdata)
      
      mapdata_validated=compute_HiC_column_in_mapdata(PCHiC, mapdata,validated,cell_type)
      mapdata$HIC=mapdata_validated
      
      crd_dist_bins = c(0,1,1e04,2e04,5e04,1e05,2e05,5e05,1e06)
      mat_hic = matrix(0,nrow=(length(crd_dist_bins)-1),ncol=1)
      for(i in 1:(length(crd_dist_bins)-1)){
        mat_hic[i,1] = compute_ratio_hic_mapdata(crd_dist_bins[i],crd_dist_bins[i+1])
      }
      
      ############# PLOT 4.3 #################

      toplot = data.frame(counts = mat_hic[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
      toplot$dist <- factor(toplot$dist ,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
      
      pdf(paste0(path_out,'/4.3_',name_condition,".pdf"))
      
      ga <- gplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene-CRD associations") +
         geom_bar(stat = "identity",fill="#56B4E9") +
         labs(x = "Distance",y = "Fraction with HiC support (%)") +
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))
      
      save(toplot,file=paste0(path_out,'4.3_',name_condition,".txt"))
      dev.off()
      
      #######################
      #Co-expression analysis
      #######################
      
      #Annotate corr_genes with CRD mapping
      
      CRD_IDs = unique(mapdata$CRD_ID)
      IDXs = c()
      CRDdistance = c()
      CRDdistance_gene1=c()
      CRDdistance_gene2=c()
      hicsupport = c()
      corr_genes$dist1=-1
      corr_genes$dist2=-1
      
      for(i in 1:length(CRD_IDs)){
        # cat(i,"\n")
        idx = which(mapdata$CRD_ID == CRD_IDs[i])
        if(length(idx)>1){
          for(k in 1:(length(idx)-1)){
          	for(l in (k+1):length(idx)){
            	geneA = mapdata$phenotype_ID[idx[k]]
            	geneB = mapdata$phenotype_ID[idx[l]]
              myhit = c(unlist(corr_genes[gene1 == geneB & gene2 == geneA,9]),unlist(corr_genes[gene1 == geneA & gene2 == geneB,9]))
              if(length(myhit)>0){
                      #cat(CRD_IDs[i],geneA,geneB,sep='\n')
                IDXs = c(IDXs,myhit)
                tmp = mean(abs(mapdata$distance[idx[k]]),abs(mapdata$distance[idx[l]])) # distance of the gene to the CRD CRD_IDs[i]
                tmp2 = mean(mapdata_validated[idx[k]],mapdata_validated[idx[l]])
                hicsupport = c(hicsupport,tmp2)
                
                CRDdistance = c(CRDdistance,tmp) # mean distance with the 2 genes
                CRDdistance_gene1 = c(CRDdistance_gene1,abs(mapdata$distance[idx[k]])) # dist of gene 1
                CRDdistance_gene2 = c(CRDdistance_gene2,abs(mapdata$distance[idx[l]])) # dist of gene 2
              }
            }
          }
        }
      }
      
      tmpdf = data.table(idx=IDXs,dist=CRDdistance,hicsupport=hicsupport, dist1=CRDdistance_gene1, dist2=CRDdistance_gene2)
      tmpdf <- tmpdf %>%
          group_by(idx) %>%
          summarize(mindist=min(dist),hic=mean(hicsupport), mindist1=min(dist1), mindist2=min(dist2))
      
      corr_genes$sameCRD[tmpdf$idx] = tmpdf$mindist # by default corr_genes$sameCRD=-1. replace it with the min distance to the CRDs
      corr_genes$hic = 0
      corr_genes$hic[tmpdf$idx] = tmpdf$hic
      corr_genes$dist1[tmpdf$idx] = tmpdf$mindist1
      corr_genes$dist2[tmpdf$idx] = tmpdf$mindist2
      
      compute_ratio <-function(coexpressed=T,mindist,maxdist,CRDmindist,CRDmaxdist,corr_genes,pval.cutoff=0.01){
          if(coexpressed){
              up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
              down = sum(corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          } else {
              up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
              down = sum(corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          }
          up/down*100
      }
      
      # new 
      compute_ratio_n <-function(coexpressed=T,mindist,maxdist,CRDmindist,CRDmaxdist,corr_genes,pval.cutoff=0.01){
        if(coexpressed){
            # corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist
            up_both = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 ==0) & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
            up_one = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 !=0) & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
            down = sum(corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          } else {
          up_both = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 ==0) & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          up_one = sum( corr_genes$sameCRD==0 & (corr_genes$dist1 + corr_genes$dist2 !=0) & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          down = sum(corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
        }
        c(up_both/down*100, up_one/down*100)
      }
      
      compute_OR <-function(mindist,maxdist,pval.cutoff=0.01){
          coexpressed_same_CRD = sum(corr_genes$sameCRD>=0 & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          coexpressed_notsame_CRD = sum(!(corr_genes$sameCRD>=0) & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          notcoexpressed_same_CRD = sum(corr_genes$sameCRD>=0 & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          notcoexpressed_notsame_CRD = sum(!(corr_genes$sameCRD>=0) & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          m = matrix(c(coexpressed_same_CRD,coexpressed_notsame_CRD,notcoexpressed_same_CRD,notcoexpressed_notsame_CRD),ncol=2)
          pval = fisher.test(m)$p.value
          OR = fisher.test(m)$estimate
          list(OR=OR,pval=pval)
      }
      
      OR_less_10kb = compute_OR(0,1e04)
      OR_less_20kb = compute_OR(1e04,2e04)
      OR_less_50kb = compute_OR(2e04,5e04)
      OR_less_100kb = compute_OR(5e04,1e05)
      OR_less_200kb = compute_OR(1e05,2e05)
      OR_less_500kb = compute_OR(2e05,5e05)
      OR_less_1Mb = compute_OR(5e05,1e06)
      OR_all = c(OR_less_10kb[[1]],OR_less_20kb[[1]],OR_less_50kb[[1]],OR_less_100kb[[1]],OR_less_200kb[[1]],OR_less_500kb[[1]],OR_less_1Mb[[1]])
      pval_all = c(OR_less_10kb[[2]],OR_less_20kb[[2]],OR_less_50kb[[2]],OR_less_100kb[[2]],OR_less_200kb[[2]],OR_less_500kb[[2]],OR_less_1Mb[[2]])
      
      
      gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06) # dist btw genes
      crd_dist_bins = c(0,1,1e04,1e05,1e06) # dist to CRD
      coexpressed_mat = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)-1))
      notcoexpressed_mat = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)-1))
      for(i in 1:(length(gene_dist_bins)-1)){
          for(k in 1:(length(crd_dist_bins)-1)){
              coexpressed_mat[i,k] = compute_ratio(T,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
              notcoexpressed_mat[i,k] = compute_ratio(F,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
          }
      }
      
      # new version
      gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06) # dist btw genes
      crd_dist_bins = c(0,1,1e04,1e05,1e06) # dist to CRD, here to change
      coexpressed_mat_n = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)))
      notcoexpressed_mat_n = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)))
      for(i in 1:(length(gene_dist_bins)-1)){
        for(k in 1:(length(crd_dist_bins)-1)){
          if (k>1){
            km=k+1
            coexpressed_mat_n[i,km] = compute_ratio(T,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1],corr_genes)   
            notcoexpressed_mat_n[i,km] = compute_ratio(F,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1],corr_genes)
          } else {
            coexpressed_mat_n[i,c(k,k+1)]=compute_ratio_n(T,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1],corr_genes)
            notcoexpressed_mat_n[i,c(k,k+1)] = compute_ratio_n(F,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1],corr_genes)
          }
        }
      }
      coexpressed_mat_n2=coexpressed_mat_n
      notcoexpressed_mat_n2=notcoexpressed_mat_n
      coexpressed_mat2=coexpressed_mat
      notcoexpressed_mat2=notcoexpressed_mat
      
      coexpressed_mat = data.frame(coexpressed_mat)
      colnames(coexpressed_mat) = c("inside","<10kb","<100kb","<1Mb")
      coexpressed_mat$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
      coexpressed_mat$row <- factor(coexpressed_mat$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
      coexpressed_table <- melt(coexpressed_mat, id.vars = "row")
      colnames(coexpressed_table) = c("Gene","CRD","value")
      
      notcoexpressed_mat = data.frame(notcoexpressed_mat)
      colnames(notcoexpressed_mat) = c("inside","<10kb","<100kb","<1Mb")
      notcoexpressed_mat$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
      notcoexpressed_mat$row <- factor(coexpressed_mat$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
      notcoexpressed_table <- melt(notcoexpressed_mat, id.vars = "row")
      colnames(notcoexpressed_table) = c("Gene","CRD","value")
      
      coexpressed_mat_n = data.frame(coexpressed_mat_n)
      colnames(coexpressed_mat_n) = c("inside2","inside1","<10kb","<100kb","<1Mb")
      coexpressed_mat_n$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
      coexpressed_mat_n$row <- factor(coexpressed_mat_n$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
      coexpressed_table_n <- melt(coexpressed_mat_n, id.vars = "row")
      colnames(coexpressed_table_n) = c("Gene","CRD","value")
      
      notcoexpressed_mat_n = data.frame(notcoexpressed_mat_n)
      colnames(notcoexpressed_mat_n) = c("inside2","inside1","<10kb","<100kb","<1Mb")
      notcoexpressed_mat_n$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
      notcoexpressed_mat_n$row <- factor(notcoexpressed_mat_n$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
      notcoexpressed_table_n <- melt(notcoexpressed_mat_n, id.vars = "row")
      colnames(notcoexpressed_table_n) = c("Gene","CRD","value")
      
      #HiC data
      
      compute_ratio_hic <-function(coexpressed=T,CRDmindist,CRDmaxdist,mindist,maxdist,pval.cutoff=0.01,threshold=5){
        if(coexpressed){
          up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$hic>threshold & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          down = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj<pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
        } else {
          up = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$hic>threshold & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
          down = sum(corr_genes$sameCRD>=CRDmindist & corr_genes$sameCRD<CRDmaxdist & corr_genes$pval.adj>=pval.cutoff & abs(corr_genes$distance)>=mindist & abs(corr_genes$distance)<maxdist,na.rm=T)
        }
        up/down*100
      }
      
      
      ################## PLOTS
      
      Paper=c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07","#972D15")
      color_palette= Paper
      # size 572 666
      pdf(paste0(path_out,'/',name_condition,"_3.7_Co-expression_analysis_CRD_gene.pdf"))
      
      ggplot(coexpressed_table_n, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Coexpressed genes") + geom_bar(stat = "identity") + 
        scale_fill_manual(values = color_palette)  + ylim(0, 65) + 
        geom_text(  aes(   y=rep (  apply(coexpressed_mat_n[,1:5],1,sum  ),5)  ,  label = c(paste0("(",round(OR_all),")"),rep("",28))  ),vjust = -.5,  size=5) + labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
      
      ggplot(notcoexpressed_table_n, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Not Coexpressed genes") + ylim(0, 10) +
        geom_bar(stat = "identity") +scale_fill_manual(values = color_palette) +
        labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
      
      dev.off()
      
      ########## Hi-C
      #HiC support by gene to gene distance
      
      gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06,3e09)
      coexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)
      notcoexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)
      
      for(i in 1:(length(gene_dist_bins)-1)){
              coexpressed_mat_hic[i,1] = compute_ratio_hic(T,0,1e09,gene_dist_bins[i],gene_dist_bins[i+1])
              coexpressed_mat_hic[i,2] = compute_ratio_hic(T,-1,0,gene_dist_bins[i],gene_dist_bins[i+1])
      }
      
      pdf(paste0(path_out,'/',name_condition,"_HiC_support_gene_pairs_with_same_CRD.pdf"))
      toplot = data.frame(counts = coexpressed_mat_hic[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
      toplot$dist <- factor(toplot$dist ,levels = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
      
      ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("_HiC support for gene pairs associated with CRD") +
         geom_bar(stat = "identity",fill="#56B4E9") +
         labs(x = "Distance between genes",y = "Fraction with HiC support (%)") +
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))
      
      dev.off()
      
      # without > 1Mb
      pdf(paste0(path_out,'/',name_condition,"_HiC_support_gene_pairs_with_same_CRD.pdf"))
      toplot = data.frame(counts = coexpressed_mat_hic[,1,-1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
      toplot=toplot[-8,]
      toplot$dist <- factor(toplot$dist ,levels = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb"))
      
      ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("_HiC support for gene pairs associated with CRD") +
        geom_bar(stat = "identity",fill="#56B4E9") +
        labs(x = "Distance between genes",y = "Fraction with HiC support (%)") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))
      
      dev.off()
      
      
      
      #HiC support by gene to CRD distance
      
      CRD_dist_bins = c(0,1,1e04,2e04,5e04,1e05,2e05,5e05,1e06)
      coexpressed_mat_hic_crddist = matrix(0,nrow=(length(CRD_dist_bins)-1),ncol=1)
      
      for(i in 1:(length(CRD_dist_bins)-1)){
              coexpressed_mat_hic_crddist[i,1] = compute_ratio_hic(T,CRD_dist_bins[i],CRD_dist_bins[i+1],0,1e09)
      }
      
      pdf(paste0(path_out,'/',name_condition,"_HiC_support_gene_pairs_with_same_CRD_geneCRD_dist.pdf"))
      toplot = data.frame(counts = coexpressed_mat_hic_crddist[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb"))
      toplot$dist <- factor(toplot$dist ,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb"))
      
      ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene pairs") +
         geom_bar(stat = "identity",fill="#56B4E9") +
         labs(x = "Distance between genes",y = "Fraction with HiC support (%)") +
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))
      
      dev.off()
      # save(PCHiC,trans_data,trans_data_q001,trans_data_q001s,hic_validated,Q,file=paste0(name,"__analysis.rda"))
      
      
    }
  }
}      
      
      
