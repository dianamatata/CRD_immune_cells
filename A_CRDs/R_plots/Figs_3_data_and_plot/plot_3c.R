
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

df.molten <- read.table("/Users/dianaavalos/Desktop/reviews_avalos/CRD_immune_cells/A_CRDs/R_plots/Figs_3_data_and_plot/3c_data.txt", sep = "\t",header=T)
Wilcoxsignif_ALL2 <- read.table("/Users/dianaavalos/Desktop/reviews_avalos/CRD_immune_cells/A_CRDs/R_plots/Figs_3_data_and_plot/3c_Wilcoxsignif_ALL_coexpressedgenes.txt", sep = "\t",header=T)
dist_subset=c("50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb")
df.molten <- df.molten %>%  filter(dist %in% dist_subset)
df.molten$value[is.na(df.molten$value)] <- 0
df.molten$dist = factor(df.molten$dist,levels = c("50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))

# Wilcoxsignif_ALL2 <- Wilcoxsignif_ALL %>%  filter(dist %in% dist_subset)
df3 = reshape2::dcast(df.molten,dist~variable)

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
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#00AFBB","#00AFBB")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "gene-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01)) +
  geom_text(data = hCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+3,
              label = format(pvalue, scientific = T,digits = 2)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# changing ylim
c <- ggplot(mCRDmolten, aes(x = dist, y = value,fill=factor(signif), alpha=abs(1-signif)))+ 
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =c("#E7B800","#E7B800"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "gene-gene distance",y = "Fraction with PCHiC support (%)") + scale_alpha(range = c(0.5, 01))+ scale_y_reverse(limits=c(45,0)) +
  geom_text(data = mCRDmoltenS,
            mapping = aes(
              x = dist,
              y = value+3,
              label = format(pvalue, scientific = T,digits = 2)))


plot_grid(a, c, ncol = 1, nrow = 2)










