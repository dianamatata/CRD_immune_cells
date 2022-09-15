
# external inputs  ---------------------------------
args = commandArgs(trailingOnly=TRUE)
COR_c1a <- loadRData(args[1])
COR_c2a <- loadRData(args[2])
N1a <- as.integer(args[3])
N2a <- as.integer(args[4])
filename=args[5]

# Packages ---------------------------------------------
library(cocor)
library(dplyr)
library(data.table) # call fread
library(R.utils) # to read .gz

# Functions  ---------------------------------

COR_c1_vs_c2 <- function(COR_c1,COR_c2,N1,N2) {
  COR_c1_vs_c2_df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(COR_c1_vs_c2_df) <- c("peak1", "peak2", "corr_c1","corr_c2", "fisher_statistic", "pvalue")
  for (i in seq(1,dim(COR_c1)[1])){
    a <- cocor.indep.groups(COR_c1$V7[i], COR_c2$V7[i], N1, N2, alternative = "two.sided",
                            test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                            data.name = NULL, var.labels = NULL, return.htest = FALSE)
    COR_c1_vs_c2_df <- rbind(COR_c1_vs_c2_df,c(COR_c1$V1[i],COR_c1$V2[i],COR_c1$V7[i],COR_c2$V7[i],a@fisher1925$statistic,a@fisher1925$p.value))
  }
  colnames(COR_c1_vs_c2_df) <- c("peak1", "peak2", "corr_c1","corr_c2", "fisher_statistic", "pvalue")
  COR_c1_vs_c2_df
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Main  ---------------------------------

compared_corr_out<- COR_c1_vs_c2(COR_c1a, COR_c2a, N1a, N2a)
print(paste0("finished ",filename))

write.table(compared_corr_out, file =filename, sep = "\t",
            row.names = TRUE, col.names = NA)


