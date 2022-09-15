## goal: stacked bar plot 3D shared genes

# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES and PATHS
#
#############################################################################################
library(UpSetR)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggmosaic)
Paper=c("#00AFBB", "#E7B800", "#FC4E07","#972D15")
color_palette= Paper

# basically grep gene column then comupte interest
directory='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/merged_TH/'
out_directory='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots/'


#############################################################################################
#
# MAIN
#
#############################################################################################

data_types = list('hist','methyl')

data_type='hist'
condition='mean'

for(data_type in data_types){

      name=paste0(data_type,'_',condition)
      cat(name)

      NEU=fread(paste0(directory,data_type,'_','neut','_',condition,'_conditional.txt.gz'))$V1
      MON=fread(paste0(directory,data_type,'_','mono','_',condition,'_conditional.txt.gz'))$V1
      TCL=fread(paste0(directory,data_type,'_','tcell','_',condition,'_conditional.txt.gz'))$V1

      N_M = length(intersect(NEU,MON))
      N_T = length(intersect(NEU,TCL))
      M_T = length(intersect(TCL,MON))
      N_M_T = length(intersect(NEU,intersect(TCL,MON)))
      
      NEU1 = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
      MON1 = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
      TCL1 = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T)

    toplot2 = data.frame(NEU=c(length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+N_T-N_M_T+0,N_M_T),
                         MON=c(length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+0+M_T-N_M_T,N_M_T),
                         TCL=c(length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T),0+N_T-N_M_T+M_T-N_M_T,N_M_T),Group =c("1cell","2cells","3cells"))
    
    
    toplot2.molten = melt(toplot2,id.vars="Group")
    colnames(toplot2.molten)[2] = "CellType"
    toplot2.molten$Group <- factor(toplot2.molten$Group ,levels = c("1cell","2cells","3cells"))
    
    ############################################################################################# PLOT1
    expressionInput <- c(NEU = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), MON = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), 
                                                   TCL = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T), 
                                                  `NEU&MON` = N_M-N_M_T, `NEU&TCL` = N_T-N_M_T,`MON&TCL` = M_T-N_M_T, `NEU&MON&TCL` = N_M_T)
    
    pdf(paste0(out_directory,name,"_Overlap_gene_sets.pdf"),paper='a4r')
    upset(fromExpression(expressionInput), order.by = "degree",text.scale=1.8)
    upset(fromExpression(expressionInput), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")
    dev.off()
    
    ############################################################################################# PLOT 2
    
    pdf(paste0(out_directory,name,"_Overlap_gene_sets_1.pdf"),paper='a4r')
    f35  <- ggplot(toplot2.molten, aes(fill=Group, y=value, x=CellType,label = scales::percent(round(value/100,digits=2)))) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Fraction Overlap") + 
      scale_fill_manual(values = color_palette) + 
      geom_text(aes(label = value),position = position_stack(vjust = .5),size = 3)+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    print(f35)
    dev.off()
    
    m=toplot2[,1:3]
    topoplot_normalized=m / t(replicate(nrow(m), colSums(m)))*100
    topoplot_normalized$Group=toplot2$Group
    topoplot_normalized.molten = melt(topoplot_normalized,id.vars="Group")
    colnames(topoplot_normalized.molten)[2] = "CellType"

    ############################################################################################# PLOT3
    
    
    pdf(paste0(out_directory,name,"_Overlap_gene_sets_norm_2.pdf"),paper='a4r')
    
    fnor  <- ggplot(topoplot_normalized.molten, aes(fill=Group, y=value, x=CellType,label = scales::percent(round(value/100,digits=2)))) +
      geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Fraction Overlap") + 
      scale_fill_manual(values = color_palette) + 
      geom_text(aes(label = round(value,digits=2)),position = position_stack(vjust = .5),size = 3)+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    print(fnor)
    dev.off()

    # mozaic plot tuto long https://www.kaggle.com/dhafer/mosaics-plots-using-ggmosaic
    pdf(paste0(out_directory,name,"_Mozaic.pdf"),paper='a4r')
    p1 <- ggplot(data = toplot2.molten) +
      geom_mosaic(aes(weight = value, x = product(CellType), fill = Group),na.rm=TRUE) + labs(x = "Cell Type",y = "Fraction Overlap", title=' ')   + 
      scale_fill_manual(values = color_palette)+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(text = element_text(size=20),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
    p1d<- ggplot_build(p1)$data %>% as.data.frame() %>% filter(.wt > 0)
    head(p1d)
    p1d$percentage=c(topoplot_normalized$NEU,topoplot_normalized$MON,topoplot_normalized$TCL)
    p1d$percentage=round(p1d$percentage,digits = 2)
    #c(55.92, 37.17, 6.91, 43.58, 47.07, 9.35, 45.78, 34.25, 19.98)
    p1d$percentage=paste0(p1d$percentage,"%")
    p2<-p1 + geom_text(data = p1d, aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2,  label = percentage))
    p2<-p2+xlab("CellType")+ylab("Fraction Overlap")
    print(p2)
    dev.off()

    
}

# toplot = data.frame(NEU1=c(55.92, 37.17, 6.91),MON1=c(43.58, 47.07, 9.35), TCL1=c(45.78, 34.25, 19.98),Group =c("1cell_type","2cell_types","3cell_types"))
