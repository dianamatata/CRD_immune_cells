library(data.table) # call fread

folder="/Users/dianaavalos/Programming/A_CRD_plots/0_CRDs/"

a = fread(file.path(folder, "quantif_M_hist_mono.bed.gz"),head = FALSE, stringsAsFactors = FALSE)
Data <- subset( a, select = -c(V1,V2,V3,V4,V5,V6) ) # drop columns 1 to 6
# hist of many rows
plot(density(as.numeric(Data[2,])))
for (row in range(3,dim(Data)[1])){ 
  lines(density(as.numeric(Data[row,])))
}


