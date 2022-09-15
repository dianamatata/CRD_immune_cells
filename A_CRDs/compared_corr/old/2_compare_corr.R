# Here we take 10 inputs and we append a line each time we run this file.
# is it slower if we need to open R to only write one line all the time?

args = commandArgs(trailingOnly=TRUE)
COR_c1 <- as.double(args[1])
pvalue_corr_c1 <- as.double(args[2])
COR_c2 <- as.double(args[3])
pvalue_corr_c2 <- as.double(args[4])
N1 <- as.integer(args[5])
N2 <- as.integer(args[6])
peaki <- args[7]
peakj <- args[8]
loci <- args[9]
locj <- args[10]
filename <- args[11]


library(cocor)

# print(paste0("R printing: ",COR_c1," ",pvalue_corr_c1," ",COR_c2," ", pvalue_corr_c2, " " ,filename, " " , peaki, " " ,peakj, " " , loci," ",locj)
a <- cocor.indep.groups(COR_c1,COR_c2, N1, N2, alternative = "two.sided",
                            test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
                            data.name = NULL, var.labels = NULL, return.htest = FALSE)

# print(c(a@fisher1925$statistic,a@fisher1925$p.value))
# print(c(COR_c1, pvalue_corr_c1, COR_c2, pvalue_corr_c2,  N1, N2, a@fisher1925$statistic,a@fisher1925$p.value))
line=c(peaki,peakj,loci,locj,COR_c1, pvalue_corr_c1, COR_c2, pvalue_corr_c2,  N1, N2, a@fisher1925$statistic,a@fisher1925$p.value)
write(line,file=filename,append=TRUE,sep = " ")
