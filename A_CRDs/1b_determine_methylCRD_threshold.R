library(mixtools)
args = commandArgs(trailingOnly=TRUE)
cell_type = args[1] # cell_type='CD19'
directory = args[2] # /home/users/a/avalosma/scratch/CEDARv1_Mvalues/2_methylCRDs

input_file=paste0(directory,'/',cell_type, '_dist.ALLchr.bed')
crdsize = scan(input_file)
crdsize.log = log10(crdsize-2)
mixmdl = normalmixEM(crdsize.log,mu = c(2.17,4.46))
thres = qnorm(0.95,mean=mixmdl$mu[1],sd=mixmdl$sigma[1])
print(thres*1000)
write(thres*1000, file = paste0(directory,"/threshold_methylCRD_",cell_type),append = FALSE, sep = " ")

# Rscript 1.0_determine_methylCRD_threshold.R 'EGAD00001002675' /home/users/a/avalosma/scratch/1_CRD
