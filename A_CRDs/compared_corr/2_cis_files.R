library(data.table) # call fread
library(R.utils) # to read .gz

folder="/home/users/a/avalosma/scratch/9_COMPARE_CORR_PEAKS/files"
for (chr in seq(1,22)){
  print(chr)
  file=paste0(folder,"/mono_neut_tcell.chr",chr,".txt")
  df = fread(file, header=FALSE) # has all chromosomes already
  df <- subset(df, df$V2 - df$V1 <= 250)
  write.table(df, file =file.path(folder, paste0("/mono_neut_tcell_CIS.chr",chr,".txt")), sep = "\t",
              row.names = FALSE, col.names =FALSE, quote=FALSE)
}
