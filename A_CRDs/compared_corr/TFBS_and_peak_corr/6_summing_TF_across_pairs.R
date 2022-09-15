library(data.table) # for data.frame
library(tidyverse) # for <- %>%
library(reshape)
library(RColorBrewer)


TF_folder="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
TF_folder="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
pairs=c("mono_neut","mono_tcell","neut_tcell")

pair="mono_neut"
index=1
threshold=10
for (index in seq(1,2)){
  for (threshold in c(10,50,100,150)){
    
    # sum over all pairs
    all_TF_info_pair=data.frame()
    for (pair in pairs){
      filename=paste0(TF_folder, "/ALL/TFs_d_",index,"_",pair,"_thre_",threshold,"ALL.txt")
      print(filename)
      df=fread(file=filename, header=TRUE)
      names <- colnames(df)
      colnames(df) <- c("TF","bin1","bin2","bin3","bin4","bin5")
      print(dim(df)[1])
      if(dim(all_TF_info_pair)[2]==0){
        all_TF_info_pair <- df} else{
          all_TF_info_pair <- rbind(all_TF_info_pair,df)
        }
    }
  
    # aggregate per bin
    b1 <- aggregate(bin1 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
    b2 <- aggregate(bin2 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
    b3 <- aggregate(bin3 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
    b4 <- aggregate(bin4 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
    b5 <- aggregate(bin5 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
    
    b <- Reduce(function(x,y) merge(x = x, y = y, by = "TF"), list(b1,b2,b3,b4,b5 ))
    
    b$sum <- rowSums( b[,2:6] )
    check_numbers <- mapply(sum,b[,-1])
    b$sum <-  b$sum/66 # 22chr 3 pairs
    out <- b[order(-b$sum),]
    
    write.table(out, file = paste0(TF_folder, "/ALL/TFs_more_represented",index,"_thre_",threshold,"ALL.txt"),
                sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)
    write.table(out$TF, file = paste0(TF_folder, "/ALL/TFs_more_represented_only",index,"_thre_",threshold,"ALL.txt"),
                sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)
    
    head(out)
    
    # PLOT
    
    bins= c(10E-10,10E-5,10E-4,10E-3,10E-2,10E-1)
    bins2 <- c("[1e-09,1e-04]","	[1e-04,10E-3]",	"[10E-3,10E-2]","	[0.01,0.1]","	[0.1,1]")
    bins2 <- c("9","4",	"3","2","1")
    
    X=seq(0, 4, 1)
    Y=c(0,0,0,0,0)
    maxY=max(out[,2:6])
    
    pdf(paste0(TF_folder,"/ALL/TFs_more_represented_",index,"_threshold_",threshold,".pdf"), 10, 10)
    
    plot(X, Y, type="n", col="black", xlim=c(min(X), max(X)), ylim=c(0, maxY+0.01), main="", xlab="-log10(P-values)", ylab="Freq", xaxt="n")
    axis(1, at=seq(0, 4, 1), labels= bins2)
    colnum=11
    mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(colnum)
    for (r in seq(10,dim(out)[1])) {
      points(X, out[r,2:6], col="grey", type="l",lwd=1)}
    for (r in seq(1,colnum-1) ){
      points(X, out[r,2:6], col=mycolors[r], type="l", lwd=2)}
    legend("topright", legend=out$TF[1:colnum-1], fill=mycolors[1:colnum-1], bty="n")
    
    dev.off()
    
  }
}

#plot(-log10(bins), -log10(MEAN$V19), xlab="-log10(Adjusted P-values from PC1)", ylab="-log10(Adjusted P-values from Mean)", main="C. Comparison of Adjusted P-values")


# compare TF most represented per pair

d1=fread(file=paste0(TF_folder, "/ALL/TFs_d_",index,"_","mono_neut","_thre_",threshold,"ALL.txt"), header=TRUE)
d2=fread(file=paste0(TF_folder, "/ALL/TFs_d_",index,"_","mono_tcell" ,"_thre_",threshold,"ALL.txt"), header=TRUE)
d3=fread(file=paste0(TF_folder, "/ALL/TFs_d_",index,"_","neut_tcell","_thre_",threshold,"ALL.txt"), header=TRUE)

d1$sum <- rowSums( d1[,2:6] )
d1 <- b[order(-d1$sum),]

d2$sum <- rowSums( d2[,2:6] )
d2 <- b[order(-d2$sum),]

d3$sum <- rowSums( d3[,2:6] )
d3 <- b[order(-d3$sum),]
# d1$TF[1:50]
# d2$TF[1:50]
# d3$TF[1:50]

sumTF <- c(d1$TF[1:50],d2$TF[1:50],d3$TF[1:50])
print(table(sumTF))
# they have the same most importants TFs  except 1 for threshold 10 index 1


