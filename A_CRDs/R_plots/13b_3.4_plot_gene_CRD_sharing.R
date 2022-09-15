# GOAL: Plot 3.4 from paper, tissue sharing gene CRDs
# update with other tissue sharing files


# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(ggplot2)
library(corrplot)
library(reshape2)
library(data.table)
library(ggmosaic)


#############################################################################################
#
# FUNCTION get_common_fraction
#
#############################################################################################

get_common_fraction <- function(set,refset,mindist,maxdist){
    if(mindist == 0 & maxdist == 0){
        set.filtered = set[which(set$distance==0),c(1,8)]
        refset.filtered = refset[which(refset$distance==0),c(1,8)]
    } else if(maxdist > 0) {
        set.filtered = set[which(set$distance>mindist & set$distance<=maxdist),c(1,8)]
        refset.filtered = refset[which(refset$distance>mindist & refset$distance<=maxdist),c(1,8)]
    } else {
        set.filtered = set[,c(1,8)]
        refset.filtered = refset[,c(1,8)]
    }
    myfraction = sum(duplicated(rbind(set.filtered,refset.filtered)))/nrow(refset.filtered)
    myfraction
}



get_fraction_c1_replicated_in_c2_binned <- function(map_cell1_vs_cell2,map_c_vs_c){
  fraction_neut_replicated_in_tcel_all = c(get_common_fraction(map_cell1_vs_cell2,map_c_vs_c,0,0),
                                           get_common_fraction(map_cell1_vs_cell2,map_c_vs_c,1,1e03),
                                           get_common_fraction(map_cell1_vs_cell2,map_c_vs_c,1e03,1e04),
                                           get_common_fraction(map_cell1_vs_cell2,map_c_vs_c,1e04,1e05),
                                           get_common_fraction(map_cell1_vs_cell2,map_c_vs_c,1e05,1e06))
}

#############################################################################################
#
# FUNCTION PLOTS
#
#############################################################################################


plot_fractions <-function(fractions,discovered,replicated) {
    filename = paste0("CRD_gene_association_discovered_in_",discovered,"_replicated_in_",replicated,".pdf")
    pdf(filename,paper="a4r")
    dat_assos = data.frame(Distance = c("0","1-\n1e03","1e03-\n1e04","1e04-\n1e05","1e05-\n1e06"),fraction = fractions)
    dat_assos$Distance <- factor(dat_assos$Distance,levels = c("0","1-\n1e03","1e03-\n1e04","1e04-\n1e05","1e05-\n1e06"))

    p <- ggplot(dat_assos, aes(x = Distance, y = fraction))+ ggtitle(paste0("Gene-CRD associations in ",discovered)) +
      geom_bar(stat = "identity",fill="steelblue") + ylim(0,0.5) +
      labs(x = "Distance between genes and CRDs (bp)",y = paste0("Fraction of significant associations in ",replicated)) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
    print(p)
    dev.off()
}


#############################################################################################
#
# plot 3.4 CRD-gene tissue sharing
#
#############################################################################################

plot_correlation_matrix_CRD_sharing <- function(Matrix, name, plot_directory){
  
  pdf(paste0(plot_directory,name,"_3.4_Gene_CRD_tissue_sharing.pdf"))
  colnames(Matrix) = c("Neutrophils","Monocytes","T cells")
  rownames(Matrix) = c("Neutrophils","Monocytes","T cells")
  corrplot(Matrix,is.corr=F,cl.lim = c(0, 1),p.mat = Matrix,sig.level=-1,insig = "p-value",number.cex=1.5)
  dev.off()
  
}


#############################################################################################
#
# DIRECTORIES AND FILES
#
#############################################################################################

# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
path_out="/Users/dianaavalos/Programming/A_CRD_plots/figs_7_Rfile/"
directory='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/CRD_genes_scripts/analysis_files/'
plot_directory="/Users/dianaavalos/Programming/A_CRD_plots/HiC13b/"


#############################################################################################
#
# MAIN
#
#############################################################################################

condition='mean'
name=paste0("hist",'_',condition)

# load data
# loop on all conditions
# TODO: need my own directory with my values!!!!!!
map70_vs_70 = read.table(paste0(directory,'70_vs_70_mapping_gene_CRD_',condition,'_ALL.txt'))
map70_vs_72 = read.table(paste0(directory,'70_vs_72_mapping_gene_CRD_',condition,'_ALL.txt'))
map70_vs_73 = read.table(paste0(directory,'70_vs_73_mapping_gene_CRD_',condition,'_ALL.txt'))
map72_vs_70 = read.table(paste0(directory,'72_vs_70_mapping_gene_CRD_',condition,'_ALL.txt'))
map72_vs_72 = read.table(paste0(directory,'72_vs_72_mapping_gene_CRD_',condition,'_ALL.txt'))
map72_vs_73 = read.table(paste0(directory,'72_vs_73_mapping_gene_CRD_',condition,'_ALL.txt'))
map73_vs_70 = read.table(paste0(directory,'73_vs_70_mapping_gene_CRD_',condition,'_ALL.txt'))
map73_vs_72 = read.table(paste0(directory,'73_vs_72_mapping_gene_CRD_',condition,'_ALL.txt'))
map73_vs_73 = read.table(paste0(directory,'73_vs_73_mapping_gene_CRD_',condition,'_ALL.txt'))



# load data
# loop on other conditions
# TODO: need my own directory with my values!!!!!!
directory='/Users/dianaavalos/Programming/A_CRD_plots/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/significants/'
condition='mean'

condition='loom'
data_type='hist'
data_type='methyl'



name=paste0('FDR_0.05',data_type,'__',condition)
# FDR_0.05_hist_tcell_vs_tcell_mean_mapping_CRD_gene_permuts.significant.txt
map70_vs_70 = read.table(paste0(directory,'FDR_0.05_',data_type,'_neut_vs_neut_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map70_vs_72 = read.table(paste0(directory,'FDR_0.05_',data_type,'_neut_vs_mono_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map70_vs_73 = read.table(paste0(directory,'FDR_0.05_',data_type,'_neut_vs_tcell_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map72_vs_70 = read.table(paste0(directory,'FDR_0.05_',data_type,'_mono_vs_neut_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map72_vs_72 = read.table(paste0(directory,'FDR_0.05_',data_type,'_mono_vs_mono_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map72_vs_73 = read.table(paste0(directory,'FDR_0.05_',data_type,'_mono_vs_tcell_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map73_vs_70 = read.table(paste0(directory,'FDR_0.05_',data_type,'_tcell_vs_neut_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map73_vs_72 = read.table(paste0(directory,'FDR_0.05_',data_type,'_tcell_vs_mono_',condition,'_mapping_CRD_gene_permuts.significant.txt'))
map73_vs_73 = read.table(paste0(directory,'FDR_0.05_',data_type,'_tcell_vs_tcell_',condition,'_mapping_CRD_gene_permuts.significant.txt'))


colnames(map70_vs_70) = colnames(map70_vs_72) = colnames(map73_vs_73) =colnames(map70_vs_73) = colnames(map72_vs_70) = colnames(map72_vs_72) = colnames(map72_vs_73) = colnames(map73_vs_70) = colnames(map73_vs_72) =  c("gene","gene_chr","gene_start","gene_end","gene_strand","dummy","distance","CRD",
"CRD_chr","CRD_start","CRD_end","pval","slope","rank")

# set = map72_vs_70
# refset = map70_vs_70
# mindist = -1
# maxdist = 1
dev.off()
# debug
fractions=fraction_neut_replicated_in_tcel_all
discovered="neutrophils"
replicated="T cells"


fraction_neut_replicated_in_mono = get_common_fraction(map72_vs_70,map70_vs_70,-1,-1)
fraction_neut_replicated_in_tcel = get_common_fraction(map73_vs_70,map70_vs_70,-1,-1)
fraction_mono_replicated_in_neut = get_common_fraction(map70_vs_72,map72_vs_72,-1,-1)
fraction_mono_replicated_in_tcel = get_common_fraction(map73_vs_72,map72_vs_72,-1,-1)
fraction_tcel_replicated_in_neut = get_common_fraction(map70_vs_73,map73_vs_73,-1,-1)
fraction_tcel_replicated_in_mono = get_common_fraction(map72_vs_73,map73_vs_73,-1,-1)

Matrix = matrix(c(1,fraction_neut_replicated_in_mono,fraction_neut_replicated_in_tcel,
                  fraction_mono_replicated_in_neut,1,fraction_mono_replicated_in_tcel,
                  fraction_tcel_replicated_in_neut,fraction_tcel_replicated_in_mono,1),ncol=3,byrow=T)

plot_correlation_matrix_CRD_sharing(Matrix, name, plot_directory)



#############################################################################################
#
# FORGET ABOUT: not in paper
#
#############################################################################################


# bar plot of the distance between genes and CRDs associations

fraction_neut_replicated_in_mono_all = get_common_fraction(map72_vs_70,map70_vs_70,-1,-1)
plot_fractions(fraction_neut_replicated_in_mono_all,"neutrophils","monocytes")

fraction_neut_replicated_in_tcel_all=get_fraction_c1_replicated_in_c2_binned(map73_vs_70,map70_vs_70)
plot_fractions(fraction_neut_replicated_in_tcel_all,"neutrophils","T cells")

fraction_mono_replicated_in_neut_all=get_fraction_c1_replicated_in_c2_binned(map70_vs_72,map72_vs_72)
plot_fractions(fraction_mono_replicated_in_neut_all,"monocytes","neutrophils")

fraction_mono_replicated_in_tcel_all=get_fraction_c1_replicated_in_c2_binned(map73_vs_72,map72_vs_72)
plot_fractions(fraction_mono_replicated_in_tcel_all,"monocytes","T cells")

fraction_tcel_replicated_in_neut_all=get_fraction_c1_replicated_in_c2_binned(map70_vs_73,map73_vs_73)
plot_fractions(fraction_tcel_replicated_in_neut_all,"T cells","neutrophils")

fraction_tcel_replicated_in_mono_all=get_fraction_c1_replicated_in_c2_binned(map72_vs_73,map73_vs_73)
plot_fractions(fraction_tcel_replicated_in_mono_all,"T cells","monocytes")

