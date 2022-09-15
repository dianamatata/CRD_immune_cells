# plots 3.5 3.6 3.7 4.1 4.3 4.4 6.5
# mozaic plots: https://edwinth.github.io/blog/ggmm/
library(tidyverse)
library(ggplot2)
library(reshape2)
library(data.table)
library(ggmosaic)

# https://github.com/karthik/wesanderson and https://wesandersonpalettes.tumblr.com/ for color palettes
BottleRocket2 = c("#FAD510", "#CB2314", "#273046")
Rushmore1 = c("#E1BD6D", "#0B775E", "#35274A" ,"#F2300F")
Royal1 = c("#899DA4", "#C93312","#74A089","#DC863B")
Zissou1 = c("#3B9AB2", "#EBCC2A", "#B40F20")
Zissou_b_j = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00")
Darjeeling1 = c("#B40F20", "#00A08A", "#F2AD00")
Moonrise2 = c("#798E87", "#C27D38", "#CCC591", "#29211F")
Cavalcanti1 = c("#D8B70A", "#81A88D", "#972D15","#02401B") #
GrandBudapest1 = c("#F1BB7B", "#D67236", "#5B1A18")
IsleofDogs2 = c("#AA9486", "#B6854D", "#39312F")
IsleofDogs2b = c("#EAD3BF", "#AA9486", "#B6854D")
Other=c("#9A8822","#DC863B","#74A089")
Paper=c("#00AFBB", "#E7B800", "#FC4E07","#972D15")

color_palette= Paper
dev.off()


### Fig 3.5
# how did i find this data? >fraction_of_genes_shared_btw_tissues.py

toplot = data.frame(NEU=c(55.92, 37.17, 6.91),MON=c(43.58, 47.07, 9.35), TCL=c(45.78, 34.25, 19.98),Group =c("1cell_type","2cell_types","3cell_types"))

toplot.molten = melt(toplot,id.vars="Group")
colnames(toplot.molten)[2] = "CellType"

#pdf("HiC_validation_ALL_CellTypes.pdf",7,5)
toplot.molten$Group <- factor(toplot.molten$Group ,levels = c("1cell_type","2cell_types","3cell_types"))
f35  <- ggplot(toplot.molten, aes(fill=Group, y=value, x=CellType,label = scales::percent(round(value/100,digits=2)))) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Trans eGenes") + 
  scale_fill_manual(values = color_palette) + 
  geom_text(aes(label = value),position = position_stack(vjust = .5),size = 3)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
print(f35)

# not normalized
toplot_not_norm = data.frame(NEU=c(2622, 1743, 324),MON=c(1571, 1697, 337), TCL=c(786, 588, 343),Group =c("1cell_type","2cell_types","3cell_types"))
toplot.molten = melt(toplot_not_norm,id.vars="Group")
colnames(toplot.molten)[2] = "CellType"

#pdf("HiC_validation_ALL_CellTypes_nn.pdf",7,5)
toplot.molten$Group <- factor(toplot.molten$Group ,levels = c("1cell_type","2cell_types","3cell_types"))
f35  <- ggplot(toplot.molten, aes(fill=Group, y=value, x=CellType,label = scales::percent(round(value/100,digits=2)))) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Trans eGenes") + 
  scale_fill_manual(values = color_palette) + 
  geom_text(aes(label = value),position = position_stack(vjust = .5),size = 3)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
print(f35)

# mozaic plot tuto long https://www.kaggle.com/dhafer/mosaics-plots-using-ggmosaic
p1 <- ggplot(data = toplot.molten) +
  geom_mosaic(aes(weight = value, x = product(CellType), fill = Group),na.rm=TRUE) + labs(x = "Cell Type",y = "Trans eGenes", title=' ')   + scale_fill_manual(values = color_palette)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
p1
p1d<- ggplot_build(p1)$data %>% as.data.frame() %>% filter(.wt > 0)
head(p1d)
p1d$numers=c(2622, 1743, 324,1571, 1697, 337,786, 588, 343)
p1d$percentage=c(55.92, 37.17, 6.91, 43.58, 47.07, 9.35, 45.78, 34.25, 19.98)
p1d$percentage=paste0(p1d$percentage,"%")
p2<-p1 + geom_text(data = p1d, aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2,  label = percentage))
p2<-p2+xlab("CellType")+ylab("Trans eGenes")
p2

# mozaic plot working, need p1d 
toplot.molten$Group <- factor(toplot.molten$Group ,levels = c("1cell_type","2cell_types","3cell_types"))

basic <- ggplot(data = toplot.molten) +
  geom_mosaic(aes(weight = value, x = product(CellType), fill = Group),na.rm=TRUE) + labs(x = "Cell Type",y = "Trans eGenes", title=' ')   + scale_fill_manual(values = color_palette)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))+
  geom_text(data = p1d, aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2,  label = percentage))



### better code avaliable in 13th Fig 3.6:  Number of genes or CRDs as a function of the number of CRDs and vice versa, respectively. - in neutrophils
# question: could we merge them in one plot?
# pdf("Connectivity_CRD_gene.pdf")
CRDhist_70_counts <-c(1555,665,305,119,139)
genehist_70counts <-c(2989,798,171, 7)
nb_CRD_not_associated_70=4883
nb_genes_not_associated_70=8732
genehist = hist(table(mapdata$phenotype_ID),breaks=c(0,1,2,3,100),plot=F)
#barplot(c(nb_CRD_not_associated_70,CRDhist_70_counts),names=c("0","1","2","3","4","5+"),main="Number of genes associated with each CRD",cex.names=1.5,cex.axis=1.5)
#barplot(c(nb_genes_not_associated_70,genehist_70counts),names=c("0","1","2","3","4+"),main="Number of CRDs associated with each gene",cex.names=1.5,cex.axis=1.5)

ggplot(data.frame(counts = c(nb_CRD_not_associated_70,CRDhist_70_counts),Number = c("0","1","2","3","4","5+")), aes(x = Number, y = counts))+ ggtitle("Genes associated with each CRD") +
  geom_bar(stat = "identity",fill=color_palette[1]) +
  geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated genes",y = "CRD counts") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))

ggplot(data.frame(counts = c(nb_genes_not_associated_70,genehist_70counts),Number = c("0","1","2","3","4+")), aes(x = Number, y = counts))+ ggtitle("CRDs associated with each gene") +
  geom_bar(stat = "identity",fill=color_palette[2]) +
  geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated CRDs",y = "Gene counts") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))



### Fig 3.7:  Percentages of pairs of co-expressed genes (5%FDR) that associate with the same CRD as a function of distance between genes.-in neutrophils

pdf("Co-expression_analysis_CRD_gene.pdf",paper="a4r")

coexpressed_table_70 <- readRDS("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/coexpressed_table_70.rds")
notcoexpressed_table_70 <- readRDS("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/notcoexpressed_table_70.rds")

ggplot(coexpressed_table_70, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Coexpressed genes") + ylim(0, 65) +
  geom_bar(stat = "identity") +
  geom_text(aes(y=rep(apply(coexpressed_mat[,1:4],1,sum),4),label = c(paste0("(",round(OR_all),")"),rep("",21))),vjust = -.5,  size=3) + labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") + scale_fill_manual(values =  color_palette) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))

ggplot(notcoexpressed_table, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Not Coexpressed genes") + ylim(0, 10) +
  geom_bar(stat = "identity") +
  labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") + scale_fill_manual(values =  color_palette) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))

dev.off()


### Fig 4.1:  Percentages of significant association at 1% FDR between pairs of peaks located on the same chromosome as a function of PCHiC signal (CHiCAGO score)
color_palette= Paper

toplot = data.frame(NEU=c(2.798961,4.452395, 5.258975, 5.542018,6.381534,6.789875,9.078163,8.726287,9.677419, 11.404562, 10.830173),
                     MON=c(2.182194,3.263072,3.906982,4.733846,5.738780,5.515588,5.747478,6.743257,7.581804,7.944915,8.346334),
                     TCL=c(1.818670,3.011890,4.143895,4.972622,5.799128,6.293706,6.786201,6.981254,7.950000,6.767411,9.672619),
                     Number = c(seq(2,20,2),">20"))

toplot.molten = melt(toplot,id.vars="Number")
toplot.molten$Number = factor(toplot.molten$Number,levels = toplot$Number)
colnames(toplot.molten)[2] = "CellType"

#pdf("HiC_validation_ALL_CellTypes.pdf",7,5)
f4.1  <- ggplot(toplot.molten, aes(x = Number, y = value,fill=CellType))+ ggtitle("HiC support for gene-CRD associations") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  color_palette) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "PCHiC Score",y = "Fraction correlated peaks (%)") 
print(f4.1)


#### Fig 4.3:  Fraction of CRD-gene associations supported by PCHiC data (CHiCAGO score >5) at increasing CRD-gene distances. 
#692, 565 good size
toplot2 = data.frame(NEU=c(33.01731, 26.62890, 24.71910, 21.97393, 23.14815, 22.64875, 15.50095, 11.81319),
                     MON=c(38.554217, 33.132530, 24.719101, 29.946524, 37.068966, 39.062500, 25.619835,  8.108108),
                     TCL=c(33.526235, 23.300971, 23.102310, 26.016260, 23.843416, 25.641026, 19.238901,  9.234828),
                     dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))


toplot2.molten = melt(toplot2,id.vars="dist")
toplot2.molten$dist = factor(toplot2.molten$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
colnames(toplot2.molten)[2] = "CellType"

f4.3  <- ggplot(toplot2.molten, aes(x = dist, y = value,fill=CellType))+ ggtitle("HiC support for gene-CRD associations") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  color_palette) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Distance",y = "Fraction with HiC support (%)") 
print(f4.3)



### Fig 4.4 Fraction of CRD-gene associations supported by PCHiC data (mean CHiCAGO score >5) for pairs of co-expressed genes (5%FDR) that associate with the same CRD


toplot2 = data.frame(NEU=c(35.95238, 38.22394, 36.69468, 43.73297, 46.57360, 43.96843, 63.79310, 77.17842),
                     MON=c(29.59184, 31.05023, 38.69732, 35.12476, 39.72366, 50.15480, 48.20847, 49.00662),
                     TCL=c(30.39216, 28.00000, 24.82759, 21.48148, 31.85841, 47.48201, 49.41176, 70.68966),
                     dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))

toplot2.molten = melt(toplot2,id.vars="dist")
toplot2.molten$dist = factor(toplot2.molten$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
colnames(toplot2.molten)[2] = "CellType"

f4.4 <- ggplot(toplot2.molten, aes(x = dist, y = value,fill=CellType))+ ggtitle("HiC support for gene pairs associated with CRD") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  color_palette) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Distance between genes",y = "Fraction with HiC support (%)") 
print(f4.4)

### Fig 6.5: Overlap between trans eGenes (aCRD, sCRD and eGene) across the three immune cell types.

toplot2 = data.frame(NEU=c(67,1,1),MON=c(39,1,1),TCL=c(3,0,1), Group =c("1cell_type","2cell_types","3cell_types"))
toplot2.molten = melt(toplot2,id.vars="Group")
colnames(toplot2.molten)[2] = "CellType"

toplot2.molten$Group <- factor(toplot2.molten$Group ,levels = c("1cell_type","2cell_types","3cell_types"))
f6.5 <-ggplot(toplot2.molten, aes(fill=Group, y=value, x=CellType,label = scales::percent(round(value/100,digits=2)))) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Trans eGenes") + 
  scale_fill_manual(values =  color_palette) +  geom_text(aes(label = round(value/100,digits=2)), position = position_stack(vjust = .5),size = 3)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
# scale_fill_brewer(palette = "Set3") +
print(f6.5)

## plots
#optional to normalize
topoplot_normalized=apply(toplot2[1:3],2,norm<-function(x){return (x/sum(x)*100)})
topoplot_normalized <- as.data.frame((topoplot_normalized))
topoplot_normalized$Group<-c("1cell_type","2cell_types","3cell_types")

toplot3.molten = melt(topoplot_normalized,id.vars="Group")
toplot3.molten$Group <- factor(toplot3.molten$Group ,levels = c("1cell_type","2cell_types","3cell_types"))
colnames(toplot3.molten)[2] = "CellType"

f65n  <-ggplot(toplot3.molten, aes(fill=Group, y=value, x=CellType, label = scales::percent(round(value/100,digits=2)))) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Trans eGenes") + 
  scale_fill_manual(values = color_palette) +  geom_text(aes(label = round(value/100,digits=2)), position = position_stack(vjust = .5),size = 3)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
print(f65n)


### mozaic plots
# https://github.com/tidyverse/ggplot2/issues/2746  https://mran.microsoft.com/snapshot/2017-08-06/web/packages/ggmosaic/vignettes/ggmosaic.html
# ggmm(toplot2.molten, CellType, Group, add_text = "prop") #blook deeper into how to use it: https://edwinth.github.io/blog/ggmm/

p <- ggplot(data = toplot2.molten) +
  geom_mosaic(aes(weight = value, x = product(CellType), fill = Group),na.rm=TRUE) + labs(x = "Cell Type",y = "Trans eGenes", title=' ')   + scale_fill_manual(values = color_palette)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
p
pd<- ggplot_build(p)$data %>% as.data.frame()
pd$percentage=paste0(round(toplot3.molten$value,digits=2),"%")
p2<-p + geom_text(data = pd, aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2,  label = percentage))
p2<-p2+xlab("CellType")+ylab("Trans eGenes")
p2


