library(data.table) # for data.frame
library(tidyverse) # for <- %>%
library(reshape)


TF_folder="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
TF_folder="/srv/beegfs/scratch/users/a/avalosma/9_COMPARE_CORR_PEAKS/TFBS_pvalue"
pairs=c("mono_neut","mono_tcell","neut_tcell")

pair="mono_neut"
index=1
threshold=10
for (pair in pairs){
  for (index in seq(1,2)){
    for (threshold in c(10,50,100,150)){
      
      all_TF_info_pair=data.frame()
      for (chr in seq(1,22)){
        filename=paste0(TF_folder,"/TFs_d_",index,"_",pair,"_",chr,"_thre_",threshold,".txt")
        if (file.exists(file=filename)){
          print(filename)
          df=fread(file=filename, header=TRUE)
          df$chr=chr
          names <- colnames(df)
          colnames(df) <- c("TF","bin1","bin2","bin3","bin4","bin5","chr")
          print(dim(df)[1])
          if(dim(all_TF_info_pair)[2]==0){
            all_TF_info_pair <- df} else{
              all_TF_info_pair <- rbind(all_TF_info_pair,df)
            }
          }
        }
        
        # aggregate per bin
        b1 <- aggregate(bin1 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
        b2 <- aggregate(bin2 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
        b3 <- aggregate(bin3 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
        b4 <- aggregate(bin4 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
        b5 <- aggregate(bin5 ~ TF, data = all_TF_info_pair, FUN = sum, na.rm = TRUE)
        
        b <- Reduce(function(x,y) merge(x = x, y = y, by = "TF"), list(b1,b2,b3,b4,b5 ))
#        colnames(b) <- colnames(df)[-7]
#        colnames(b) <- names[-7]
        
        write.table(b,
          file = paste0(TF_folder, "/TFs_d_",index,"_",pair,"_thre_",threshold,"ALL.txt"),
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE,
          quote = FALSE)
    }
  }
}

