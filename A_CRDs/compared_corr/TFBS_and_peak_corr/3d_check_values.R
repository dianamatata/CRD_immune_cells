library(data.table) # for data.frame
library(tidyverse) # for <- %>%
library(ggplot2) 

dir="/Users/dianaavalos/Programming/A_CRD_plots/9_COMPARE_CORR_PEAKS/"
dir="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr/TFBS_and_peak_corr"
df=fread(paste0(dir,"/3c_check_sizes_intersect.txt"), head=FALSE)

df$chr<-sapply(df$V1, function(x) str_split(x,"_")[[1]][3])
df$threshold <- sapply(df$V1, function(x) strtoi(str_split(str_split(x,"_")[[1]][5],"[.]")[[1]][1]))

colnames(df) <- c("file","count","chr","threshold")
#df$threshold <- as.character(df$threshold)

head(df)

df$threshold <- factor(df$threshold ,levels = c(10,50,100,150))
df$chr <- factor(df$chr ,levels = seq(1,22))

pdf(paste0(dir,"/TFs_per_chr_per_threshold.pdf"), 20, 10)

ggplot(data = df, aes(x=chr , y = count, fill=threshold)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+ labs(title="Number of TFs per chromosome per histone peak width/2",
  y ="Number of TFs", x = "Chromosomes")

dev.off()

