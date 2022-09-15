# Packages ---------------------------------------------

library(ggplot2)
library(data.table) # for data.frame
library(tidyverse) # for <- %>% 
library(qvalue)

# Functions ---------------------------------------------

image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

# Main ---------------------------------------------


folder="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TF_arrays"
folder="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TF_arrays"

pairs=c("mono_neut","mono_tcell","neut_tcell")

pair="mono_neut"
size=50
threshold=0.05

for (pair in pairs){
  print(pair)
  df_NOT_specific_NO_overlap <- read.csv(paste0(folder, "/pairs/","ALL_NOT_spec_NO_O_", pair, "_ALL.txt"), row.names = 1, header= TRUE, sep="\t")
  df_specific_NO_overlap <- read.csv(paste0(folder, "/pairs/","ALL_spec_NO_O_", pair, "_ALL.txt"), row.names = 1, header= TRUE, sep="\t")
  df_NOT_specific <- read.csv(paste0(folder, "/pairs/","ALL_NOT_spec_", pair, "_ALL.txt"), row.names = 1, header= TRUE, sep="\t")
  df_specific <- read.csv( paste0(folder, "/pairs/","ALL_spec_", pair, "_ALL.txt"), row.names = 1, header= TRUE, sep="\t")
  
  ODDS = matrix(0, nrow=size, ncol=size)
  PVAL = matrix(0, nrow=size, ncol=size)
  ODDS_PVAL = matrix(0, nrow=size, ncol=size)
  
    for (i in seq(1,size)){
      for (j in seq(1,size)){
        FT = fisher.test(matrix(c(df_specific[i,j], df_NOT_specific[i,j], df_specific_NO_overlap[i,j], df_NOT_specific_NO_overlap[i,j]), nrow=2))
        ODDS[i,j] = FT$estimate
        PVAL[i,j] = FT$p.value
        if(FT$p.value<=threshold){
          ODDS_PVAL[i,j] = FT$estimate
        } else {
          ODDS_PVAL[i,j] = 1 # so it is displayed as white on the plot
        }
      }
    }
    
    ODDS_PVAL[ODDS_PVAL >= 2] = 2

  
  write.table(ODDS,
              file = paste0(folder, "/ODDS_", pair, ".txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  write.table(PVAL,
              file = paste0(folder, "/PVAL_of_ODDS", pair, ".txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  write.table(ODDS_PVAL,
              file = paste0(folder, "/ODDS_PVAL_", pair, ".txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE) 
  # Plot ---------------------------------------------
}


for (pair in pairs){
  print(pair)

  ODDS_PVAL <- read.csv(paste0(folder, "/ODDS_PVAL_", pair, ".txt"), row.names = 1, header= TRUE, sep="\t")
  colnames(ODDS_PVAL) <- colnames(df_NOT_specific)
  row.names(ODDS_PVAL)<- row.names(df_NOT_specific)

  # Plot ---------------------------------------------
  
  pdf(paste0(folder,"/",pair,"_ALL.pdf"), 10, 10)
  
  layout(matrix(c(1,2), nrow=2, ncol=1), widths=c(8,1), heights=c(8,1)) #layout.show(2)
  breaks <- seq(0, 2,length.out=100)
  breaks <- seq(1-delta, 1+delta,length.out=100)
  
  par(mar=c(5,5,2,2))
  image(seq(dim(ODDS_PVAL)[1]), seq(dim(ODDS_PVAL)[2]), data.matrix(ODDS_PVAL),
        col=pal.1(length(breaks)-1), breaks=breaks, xaxt="n", yaxt="n", ylab="", xlab="")
  axis(1, at = 1:ncol(ODDS_PVAL), labels=colnames(df_specific_NO_overlap),cex.axis=0.6, las=2)
  axis(2, at = 1:nrow(ODDS_PVAL), labels=colnames(df_specific_NO_overlap),cex.axis=0.6, las=2)
  
  par(mar=c(3,1,1,1))
  image.scale(ODDS_PVAL, col=pal.1(length(breaks)-1), breaks=breaks, horiz=TRUE)
  
  dev.off()
  
}
