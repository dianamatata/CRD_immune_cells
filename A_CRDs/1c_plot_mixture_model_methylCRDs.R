library(mixtools)
args = commandArgs(trailingOnly=TRUE)
cell_type = args[1] # cell_type='methyl_neut'
directory = args[2] # /home/users/a/avalosma/scratch/1_CRD
outdirectory= args[3] # same as directory? # A_CRD_plots

input_file=paste0(directory,'/',cell_type, '_dist.ALLchr.bed')
crdsize = scan(input_file)
crdsize.log = log10(crdsize-2)
mixmdl = normalmixEM(crdsize.log,mu = c(2.17,4.46))
thres = qnorm(0.95,mean=mixmdl$mu[1],sd=mixmdl$sigma[1])
print(thres*1000)
write(thres*1000, file = paste0(directory,"/threshold_methylCRD_",cell_type),append = FALSE, sep = " ")

pdf(paste0(outdirectory,"/mixture_model_",cell_type,".pdf"))
plot(mixmdl,density = TRUE, w = 2,xlab2="log10(Size)",n=40)
lines(density(crdsize.log), lty=2, lwd=2)
abline(v=thres,lty=2,lwd=3)
text(x=thres+1,y=0.3,paste0("Size threshold:\n",round((10^thres)/1000,1)," kb"))
dev.off()

# Rscript 1c_plot_mixture_model_methylCRDs.R 'methyl_neut' /home/users/a/avalosma/scratch/1_CRD /home/users/a/avalosma/scratch/A_CRD_plots
