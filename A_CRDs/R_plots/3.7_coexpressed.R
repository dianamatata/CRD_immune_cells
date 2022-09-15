# 3.7 coexpressed

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


#############################################################################################
#
# Co-expression analysis
#
#############################################################################################



data_type='methyl'
cell_type='neut'
condition='loom'
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
      array_aCRD_gene = read.table(file_mapdata, hea=F, stringsAsFactors=F)
      corr_genes=fread(paste0(path,rna_file[[cell_type]]))
      corr_genes=get_corr_genes_formated(corr_genes,genelist)
      
      colnames(array_aCRD_gene) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                            "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","rank",
                            "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")

      
      #Annotate corr_genes with CRD mapping
      
      CRD_IDs = unique(array_aCRD_gene$CRD_ID)
      IDXs = c()
      CRDdistance = c()
      CRDdistance_gene1=c()
      CRDdistance_gene2=c()
      hicsupport = c()
      corr_genes$dist1=-1
      corr_genes$dist2=-1
      
      for(i in 1:length(CRD_IDs)){
        # cat(i,"\n")
        idx = which(array_aCRD_gene$CRD_ID == CRD_IDs[i])
        if(length(idx)>1){
          for(k in 1:(length(idx)-1)){
            for(l in (k+1):length(idx)){
              geneA = array_aCRD_gene$phenotype_ID[idx[k]]
              geneB = array_aCRD_gene$phenotype_ID[idx[l]]
              myhit = c(unlist(corr_genes[gene1 == geneB & gene2 == geneA,9]),unlist(corr_genes[gene1 == geneA & gene2 == geneB,9]))
              if(length(myhit)>0){
                #cat(CRD_IDs[i],geneA,geneB,sep='\n')
                IDXs = c(IDXs,myhit)
                tmp = mean(abs(array_aCRD_gene$distance[idx[k]]),abs(array_aCRD_gene$distance[idx[l]])) # distance of the gene to the CRD CRD_IDs[i]
      
                CRDdistance = c(CRDdistance,tmp) # mean distance with the 2 genes
                CRDdistance_gene1 = c(CRDdistance_gene1,abs(array_aCRD_gene$distance[idx[k]])) # dist of gene 1
                CRDdistance_gene2 = c(CRDdistance_gene2,abs(array_aCRD_gene$distance[idx[l]])) # dist of gene 2
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
      corr_genes$dist1[tmpdf$idx] = tmpdf$mindist1
      corr_genes$dist2[tmpdf$idx] = tmpdf$mindist2
      
      #############################################################################################
      #
      # NEXT
      #
      #############################################################################################
      
      
      
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
      
      
      #############################################################################################
      #
      # PLOT
      #
      #############################################################################################
      # need coexpressed_table_n notcoexpressed_table_n OR_all
      
      # size 572 666
      pdf(paste0(path_out,'/3.7_',name_condition,"_3.7_Co-expression_analysis_CRD_gene.pdf"))
      
      ggplot(coexpressed_table_n, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Coexpressed genes") + geom_bar(stat = "identity") + 
        scale_fill_manual(values = color_palette)  + ylim(0, 65) + 
        geom_text(  aes(   y=rep (  apply(coexpressed_mat_n[,1:5],1,sum  ),5)  ,  label = c(paste0("(",round(OR_all),")"),rep("",28))  ),vjust = -.5,  size=5) + labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
      
      ggplot(notcoexpressed_table_n, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Not Coexpressed genes")  +  geom_bar(stat = "identity") +
        scale_fill_manual(values = color_palette) + ylim(0, 10) +
        labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
        geom_text(  aes(   y=rep (  apply(coexpressed_mat_n[,1:5],1,sum  ),5)  ,  label = c(paste0("(",round(OR_all),")"),rep("",28))  ),vjust = -.5,  size=5) + labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
      
      dev.off()
      
      
    }
  }
}







