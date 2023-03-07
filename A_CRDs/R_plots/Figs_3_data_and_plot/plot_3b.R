# plotCRD_gene_hic

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

path_out="/Users/dianaavalos/Desktop/reviews_avalos/CRD_immune_cells/A_CRDs/R_plots/Figs_3_data_and_plot"
toplot_ALL=read.table(file = paste0(path_out,"/toplot_ALL_CRDgene.txt"), sep = "\t",header=T)
Wilcoxsignif_ALL=read.table(file = paste0(path_out,"/Wilcoxsignif_ALL_CRDgene.txt"), sep = "\t",header=T)

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

hCRDmolten <- df.molten %>% filter(CRD_type=="hCRD")
mCRDmolten <- df.molten %>% filter(CRD_type=="mCRD")
hCRDmolten$pvalue=rep(Wilcoxsignif_ALL[Wilcoxsignif_ALL$data_type=="hist",]$pvalue,2)
mCRDmolten$pvalue=rep(Wilcoxsignif_ALL[Wilcoxsignif_ALL$data_type=="methyl",]$pvalue,2)
mCRDmoltenS <- mCRDmolten %>% filter(variable=="mCRD_s")
hCRDmoltenS <- hCRDmolten %>% filter(variable=="hCRD_s")
for ( i in seq(1,length(df3$hCRD_s))){
  hCRDmoltenS$value[i]=max(df3$hCRD_s[i],df3$hCRD_ns[i])
}


#############################################################################################
#
# PLOT with some columns
#
#############################################################################################

df.molten <- df.molten %>%  filter(dist %in% c("20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
df3 <- df3 %>%  filter(dist %in% c("20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
Wilcoxsignif_ALL2 <- Wilcoxsignif_ALL %>%  filter(dist %in% c("20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))

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


#############################################################################################
#
# PLOT with all columns
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


