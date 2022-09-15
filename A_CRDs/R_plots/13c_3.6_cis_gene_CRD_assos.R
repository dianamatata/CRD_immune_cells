#### CRD gene associations, no HIC

# hist_neut_mean_mapping_CRD_gene_permuts.txt.gz
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
library(tidyverse)

#############################################################################################
#
# FUNCTIONS
#
#############################################################################################

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

color_palette=c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07","#972D15")


data_type='hist'
cell_type='tcell'
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
      mapdata = read.table(file_mapdata, hea=F, stringsAsFactors=F)
      corr_genes=fread(paste0(path,rna_file[[cell_type]]))
      corr_genes=get_corr_genes_formated(corr_genes,genelist)
      
      colnames(mapdata) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                  "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","rank",
                                  "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
      array_aCRD_gene=mapdata
      
      #Filter protein coding and long nc RNA genes
      array_aCRD_gene_unfiltered=array_aCRD_gene
      array_aCRD_gene = array_aCRD_gene[array_aCRD_gene$phenotype_ID %in% genelist,]
      # order by chromosome and position
      array_aCRD_gene = array_aCRD_gene[order(array_aCRD_gene[,2],array_aCRD_gene[,3]),]
      nb_CRD_not_associated = nrow(allCRDs) - length(unique(array_aCRD_gene$CRD_ID))
      nb_genes_not_associated = length(unique(corr_genes$gene1)) - length(unique(array_aCRD_gene$phenotype_ID))

    ######## PLOT
      CRDhist = hist(table(array_aCRD_gene$CRD_ID),breaks=c(0,1,2,3,4,100),plot=F)
      genehist = hist(table(array_aCRD_gene$phenotype_ID),breaks=c(0,1,2,3,100),plot=F)
      
      pdf(paste0(path_out,'/',name_condition,"_3.6_Connectivity_CRD_gene.pdf"))

      ggplot(data.frame(counts = c(nb_CRD_not_associated,CRDhist$counts),Number = c("0","1","2","3","4","5+")), aes(x = Number, y = counts))+ ggtitle("Genes associated with each CRD") +
        geom_bar(stat = "identity",fill="#E7B800") +
        geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated genes",y = "CRD counts") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
      
      ggplot(data.frame(counts = c(nb_genes_not_associated,genehist$counts),Number = c("0","1","2","3","4+")), aes(x = Number, y = counts))+ ggtitle("CRDs associated with each gene") +
        geom_bar(stat = "identity",fill="#046C9A") +
        geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated CRDs",y = "Gene counts") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
      
     dev.off()
      
    }
  }
}


colnames(array_aCRD_gene) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                              "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","degree_freedom",
                              "Dummy","1st_param_beta","2nd_param_beta","nominal_pval","r_squared","slope","empirical_pval","beta_pval")
      