# get chr 5: 1000 to 2000 bp

# Clean environment ------------------------------------

rm(list = ls())
gc()

# Packages ---------------------------------------------

library("biomaRt")
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gggenes)
library(RColorBrewer)

# Load data  ---------------------------------------------


path="/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/significants"
histCRDgene_mono <- read.table(file.path(path,"FDR_0.05_hist_mono_mean_mapping_CRD_gene_ALL.significant.txt"),header = F)
histCRDgene_neut <- read.table(file.path(path,"FDR_0.01_hist_neut_mean_mapping_CRD_gene_ALL.significant.txt"), header = F)
histCRDgene_tcell <- read.table(file.path(path,"FDR_0.01_hist_tcell_mean_mapping_CRD_gene_ALL.significant.txt"), header = F)
outpath="/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots"

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Functions  ---------------------------------

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
    mart = human,
    useCache = FALSE
  )
  gene_coords$size = gene_coords$end_position - gene_coords$start_position
  gene_coords
}

get_genes_ENSG_ID <- function(array) {
  c <- strsplit(array$V1, '\\.')
  Genes <- unlist(c)[2*(1:length(c))-1]
  Genes
}
# Parameters ---------------------------------------------


chromosome=5
start_range=39420213 
stop_range=60728752
# range for fig:1000 to 2000 peaks, we want bp
# [avalosma@login2]~/scratch/8_PEAKS% cat hist_mono.corr.chr5.txt | grep '^1000 2000'
# 1000 2000 39420897 60728190 H3K4me1_chr5-49894 H3K27ac_chr5-41819 0.041 0.604177

# Main ---------------------------------------------

# filter chromosome 5 and start stop range in bp for the peak coordinates (peak 1000 and 2000)
histCRDgene_mono %>% filter(V2 ==5)
genes_mono <- histCRDgene_mono %>% filter(V2 ==5) %>% filter( (V3 >= start_range & V3 <= stop_range) | (V4 >= start_range & V4 <= stop_range) )
genes_neut <- histCRDgene_neut %>% filter(V2 ==5) %>% filter( (V3 >= start_range & V3 <= stop_range) | (V4 >= start_range & V4 <= stop_range) )
genes_tcell <- histCRDgene_tcell %>% filter(V2 ==5) %>% filter( (V3 >= start_range & V3 <= stop_range) | (V4 >= start_range & V4 <= stop_range) )

# look for these genes on human biomart, get symbols and length  ---------------------------------------------
genes_mono$ensembl_gene_id <- get_genes_ENSG_ID(genes_mono)
genes_neut$ensembl_gene_id <- get_genes_ENSG_ID(genes_neut)
genes_tcell$ensembl_gene_id <- get_genes_ENSG_ID(genes_tcell)

gene_coords_mono = get_genes_info(human, genes_mono$ensembl_gene_id)
gene_coords_neut = get_genes_info(human, genes_neut$ensembl_gene_id)
gene_coords_tcell = get_genes_info(human, genes_tcell$ensembl_gene_id)

genes_mono_full <- left_join(genes_mono, gene_coords_mono, "ensembl_gene_id")
genes_neut_full <- left_join(genes_neut, gene_coords_neut, "ensembl_gene_id")
genes_tcell_full <- left_join(genes_tcell, gene_coords_tcell, "ensembl_gene_id")

genes_mono_full2 <- genes_mono_full[c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "strand", "size","V8","V10","V11")]
colnames(genes_mono_full2) <- c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "strand", "size","CRD","CRD_start","CRD_end")
genes_mono_full2$CRD_length <-  round((genes_mono_full2$CRD_end-genes_mono_full2$CRD_start)/10000,1)
genes_mono_full2

# need gene ensembl hgsymbol chromosome 

# plot the genes on an axis as boxes with start and stop position   ---------------------------------------------
# plot genes: https://figshare.com/articles/presentation/How_to_plot_genes_using_R/12178902/1
# https://cran.r-project.org/web/packages/gggenes/readme/README.html

# Monocytes plot    ---------------------------------------------
gene_coords_mono$molecule="DNA"
gene_coords_mono_nonan  <-  gene_coords_mono %>% filter(hgnc_symbol != "")
gene_coords_mono_nonan$end_position <- gene_coords_mono_nonan$end_position+gene_coords_mono_nonan$size # longer to be visible


ggplot(
  gene_coords_mono_nonan,
  aes(xmin = start_position, xmax = end_position, y = molecule, fill = hgnc_symbol, label = hgnc_symbol)
) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(1, "mm")) +
  #geom_gene_label(align = "left") +
  geom_blank(data = gene_coords_mono) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(15))  +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes()

ggsave(paste0(outpath,"gene_coords_mono_nonan1",".tiff"), units="in", width=20, height=4, dpi=300, compression = 'lzw')

# Neut plot    ---------------------------------------------
gene_coords_neut$molecule="DNA"
gene_coords_neut_nonan  <-  gene_coords_neut %>% filter(hgnc_symbol != "")
gene_coords_neut_nonan$end_position <- gene_coords_neut_nonan$end_position+gene_coords_neut_nonan$size


ggplot(
  gene_coords_neut_nonan,
  aes(xmin = start_position, xmax = end_position, y = molecule, fill = hgnc_symbol, label = hgnc_symbol)
) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  #geom_gene_label(align = "left") +
  #geom_blank(data = gene_coords_neut) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(16))  +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes()

ggsave(paste0(outpath,"gene_coords_neut_nonan1",".tiff"), units="in", width=20, height=4, dpi=300, compression = 'lzw')


# TCL plot    ---------------------------------------------
gene_coords_tcell$molecule="DNA"
gene_coords_tcell_nonan  <-  gene_coords_tcell %>% filter(hgnc_symbol != "")
gene_coords_tcell_nonan$end_position <- gene_coords_tcell_nonan$end_position+gene_coords_tcell_nonan$size
  
ggplot(
  gene_coords_tcell_nonan,
  aes(xmin = start_position, xmax = end_position, y = molecule, fill = hgnc_symbol, label = hgnc_symbol)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  #geom_gene_label(align = "left") +
  geom_blank(data = gene_coords_tcell) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(8))  +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() 
  #geom_text(aes(label = hgnc_symbol, x=start_position, y=molecule), vjust = -2)

ggsave(paste0(outpath,"gene_coords_tcell_nonan1",".tiff"), units="in", width=20, height=4, dpi=300, compression = 'lzw')



# which gene is linked to which CRD?
