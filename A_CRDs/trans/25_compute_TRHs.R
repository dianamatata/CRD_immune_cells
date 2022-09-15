# Goal: plot 50 TRH index for 3 cell types
# so far  this one

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
library(gplots)
library(data.table)
library(tidyverse)
library(igraph)


#############################################################################################
#
# FUNCTION
#
#############################################################################################
get_TRH <- function(TRH_signif){
  
  ###### CALL TRHs using greedy algorithm
  LINKS.TRH = data.frame(from=TRH_signif$id1, to=TRH_signif$id2, weigth=ifelse(TRH_signif$corr>0, 1, 2), stringsAsFactors=FALSE)
  NODES.TRH = data.frame(id=names(table(c(LINKS.TRH$from, LINKS.TRH$to))), chr=matrix(unlist(strsplit(names(table(c(LINKS.TRH$from, LINKS.TRH$to))), split="_")), ncol=3, byrow=TRUE)[, 1], stringsAsFactors=FALSE)
  DATAnet.TRH = graph_from_data_frame(d=LINKS.TRH, vertices=NODES.TRH, directed=FALSE)
  # do not plot takes a long time and aborts session #plot(DATAnet.TRH, vertex.label = V(DATAnet.TRH)$name)
  communities  = fastgreedy.community(DATAnet.TRH) # https://www.rdocumentation.org/packages/igraph/versions/0.4.1/topics/fastgreedy.community
  # this function tries to find dense subgraph, also called communities in graphs via directly optimizing a modularity score.
  communities
}


#############################################################################################
#
# PLOTS
#
#############################################################################################

plot_50_biggest_TRH_sizes <- function(communities,file, name){
  
  N_TRH = max(communities$membership)
  N_TRH_TOSHOW = 50
  df = data.frame(ID=1:N_TRH,SIZE=rle(sort(communities$membership))$length,stringsAsFactors=F)
  
  cat(length(df$SIZE))
  df <-df[order(-df$SIZE),]
  df <- data.frame(head(df$SIZE,50),seq(1,50))
  colnames(df)=c("SIZE","ID")
  if (length(df$SIZE) >=50){
    outfilename=paste0(out_dir,"plot_50_biggest_TRHs_",name,".pdf")
    pdf(outfilename)
    p2<- ggplot(df, aes(x = ID, y = SIZE)) + geom_bar(stat = "identity", fill="steelblue", width=0.7) + 
      scale_y_continuous(trans = 'log10') +
      labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH,"  ", name), x ="TRH index", y = "TRH Size") +
      theme(axis.text=element_text(size=16), axis.title=element_text(size=16) ) + theme_minimal() 
    print(p2)
    dev.off()
  }
  if 
  (length(df$SIZE) < 50){
    x1 <- seq(32,50)                 # Column 1 of data frame
    x2 <- rep(0,19)              # Column 2 of data frame
    df2 <- data.frame(x1, x2) 
    colnames(df2)=c("SIZE","ID")
    df3=rbind(df,df2)
    df4 <- data.frame(head(df3$ID,50),seq(1,50))
    colnames(df4)=c("SIZE","ID")
    outfilename=paste0(out_dir,"plot_biggest_TRHs_",name,".pdf")
    pdf(outfilename)
    p2<- ggplot(df4, aes(x = ID, y = SIZE)) + geom_bar(stat = "identity", fill="steelblue", width=0.7) + 
      labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH,"  ", name), x ="TRH index", y = "TRH Size") +
      theme(axis.text=element_text(size=16), axis.title=element_text(size=16) ) + theme_minimal() 
    print(p2)
    dev.off()
  }
  
}


#############################################################################################
#
# DIRECTORIES AND FILES
#
#############################################################################################

directory='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant'
out_dir='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:TRHs/'

files <- list.files(path=directory, pattern="0.0*.txt", full.names=TRUE, recursive=FALSE)
filename=paste0(out_dir,'TRHs_inventory.txt')


#############################################################################################
#
# MAIN
#
#############################################################################################
file.create(filename)
line='nbr_TRHs nbr_CRDs_involved'
write(line,file=filename,append=TRUE)

files= intersect(list.files(path=directory, pattern="0.01.txt", full.names=TRUE, recursive=FALSE), 
                 list.files(path=directory, pattern="mean", full.names=TRUE, recursive=FALSE))

for (f in files){
  # f=files[3]
  file=basename(f)
  TRH_signif = as.data.frame(data.table::fread(f, head=TRUE, stringsAsFactors=FALSE))
  colnames(TRH_signif) = c("idx1","chr1","midplace","id1","idx2","chr2","midplace2","id2","corr","pval","qvalue")
  communities=get_TRH(TRH_signif)
  nbr_TRHs = max(communities$membership)  # number of TRH
  nbr_CRDs_involved=length(communities$names) # number of CRDs involved
  COM = data.frame(communities$names, communities$membership)
  COM = COM[order(communities$membership),]
  # write info
  # line=paste0(file,' ',nbr_TRHs,' ',nbr_CRDs_involved)
  # write(line,file=filename,append=TRUE)
  name=paste0(str_sub(file, 1, - 27) ,str_sub(file, 34, - 5))
  cat (' ',name, '  ', nbr_TRHs, '  ')
  
  plot_50_biggest_TRH_sizes(communities,file, name)
}



# debug:
# f='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant/hist_mono_mean_trans.significant_0.01.txt'

