#### figs (4.1 no data) 4.3 4.4
# find where the data comes from


### Fig 4.1:  Percentages of significant association at 1% FDR between pairs of peaks located on the same chromosome as a function of PCHiC signal (CHiCAGO score)
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


#############################################################################################
#
# 4.3
#
#############################################################################################

#### Fig 4.3:  Fraction of CRD-gene associations supported by PCHiC data (CHiCAGO score >5) at increasing CRD-gene distances. 

Paper=c("#00AFBB", "#E7B800","#FC4E07")
color_palette= Paper
path_out = '/Users/dianaavalos/Programming/A_CRD_plots/figs_geneCRD/'
data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
conditions = list('mean', 'loom')


for(data_type in data_types){

      condition='mean'
      name_condition=paste0(data_type,'_',condition)
      NEU=read.delim(file=paste0(path_out,'4.3_geneCRD_associations_HiC_',data_type,'_','neut' ,'_',condition,'.txt'), header = TRUE, sep = "\t", dec = ".")
      MON=read.delim(file=paste0(path_out,'4.3_geneCRD_associations_HiC_',data_type,'_','mono' ,'_',condition,'.txt'), header = TRUE, sep = "\t", dec = ".")
      TCL=read.delim(file=paste0(path_out,'4.3_geneCRD_associations_HiC_',data_type,'_','tcell' ,'_',condition,'.txt'), header = TRUE, sep = "\t", dec = ".")

  toplot2 = data.frame(NEU=NEU$counts,MON=MON$counts,TCL=TCL$counts,
                       dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
  
  toplot2.molten = melt(toplot2,id.vars="dist")
  toplot2.molten$dist = factor(toplot2.molten$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
  colnames(toplot2.molten)[2] = "CellType"
  
  pdf(paste0(path_out,'4.3_3cells_',name_condition,".pdf"))
  
  f4.3  <- ggplot(toplot2.molten, aes(x = dist, y = value,fill=CellType))+ ggtitle("HiC support for gene-CRD associations") +
    geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  color_palette) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
    theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "Distance",y = "Fraction with HiC support (%)") 
  print(f4.3)
  
  dev.off()
  
}



# toplot2 = data.frame(NEU=c(33.01731, 26.62890, 24.71910, 21.97393, 23.14815, 22.64875, 15.50095, 11.81319),
#                      MON=c(38.554217, 33.132530, 24.719101, 29.946524, 37.068966, 39.062500, 25.619835,  8.108108),
#                      TCL=c(33.526235, 23.300971, 23.102310, 26.016260, 23.843416, 25.641026, 19.238901,  9.234828),
#                      dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))




### Fig 4.4 Fraction of CRD-gene associations supported by PCHiC data (mean CHiCAGO score >5) for pairs of co-expressed genes (5%FDR) that associate with the same CRD


for(data_type in data_types){
  condition='mean'
  name_condition=paste0(data_type,'_',condition)
  NEU=read.delim(file=paste0(path_out,'4.4_geneCRD_associations_HiC_',data_type,'_','neut' ,'_',condition,'.txt'), header = TRUE, sep = "\t", dec = ".")
  MON=read.delim(file=paste0(path_out,'4.4_geneCRD_associations_HiC_',data_type,'_','mono' ,'_',condition,'.txt'), header = TRUE, sep = "\t", dec = ".")
  TCL=read.delim(file=paste0(path_out,'4.4_geneCRD_associations_HiC_',data_type,'_','tcell' ,'_',condition,'.txt'), header = TRUE, sep = "\t", dec = ".")
  
  toplot2 = data.frame(NEU=NEU$counts,MON=MON$counts,TCL=TCL$counts,
                       dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
  
  toplot2.molten = melt(toplot2,id.vars="dist")
  toplot2.molten$dist = factor(toplot2.molten$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
  colnames(toplot2.molten)[2] = "CellType"
  
  
  pdf(paste0(path_out,'4.4_3cells_',name_condition,".pdf"))
  
  f4.4 <- ggplot(toplot2.molten, aes(x = dist, y = value,fill=CellType))+ ggtitle("HiC support for gene pairs\nassociated with CRD ") +
    geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  color_palette) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
    theme(text = element_text(size=18),axis.title = element_text(size = 18),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "Distance between co-expressed genes",y = "Fraction with HiC support (%)") +
    theme(plot.title = element_text(hjust = 0.5))
  print(f4.4)
  
  dev.off()

}


# toplot2 = data.frame(NEU=c(35.95238, 38.22394, 36.69468, 43.73297, 46.57360, 43.96843, 63.79310, 77.17842),
#                      MON=c(29.59184, 31.05023, 38.69732, 35.12476, 39.72366, 50.15480, 48.20847, 49.00662),
#                      TCL=c(30.39216, 28.00000, 24.82759, 21.48148, 31.85841, 47.48201, 49.41176, 70.68966),
#                      dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
# 
