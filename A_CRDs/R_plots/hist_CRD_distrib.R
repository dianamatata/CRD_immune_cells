
#############################################################################################
#
# histogram CRDs
#
#############################################################################################


directory='/Users/dianaavalos/Programming/A_CRD_plots/1_CRD/'
out_directory='Programming/A_CRD_plots/1_CRD/plots/'

files=list.files(path = directory, pattern = ".gz", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

for (file in files){
  crd_file=as.data.frame(data.table::fread(paste0(directory,file), head=TRUE, stringsAsFactors=FALSE))
  
  crdsize=crd_file$end-crd_file$start
  median(crdsize)
  median(crdsize.log)
  crdsize.log = log10(crdsize-2)
  
  pdf(paste0(out_directory,"histogram_",file,".pdf"))
  
  hist(crdsize.log,xlab="log10(Size)",main=paste0("median: ",round(median(crdsize)/1000,1)," kb"), xlim=c(3,8), breaks=40)
  lines(density(crdsize.log), lty=2, lwd=2)
  abline(v=  median(crdsize.log) ,lty=2,lwd=3, col='red')
  # text(x=median(crdsize.log)+1,y=max(crdsize.log)*40,paste0("median:\n",round(median(crdsize)/1000,1)," kb"), col='red')
  
  dev.off()
}
  
#############################################################################################
#
# histogram methyl CRDs
#
#############################################################################################

directory='/Users/dianaavalos/Programming/A_CRD_plots/0_CRDs/'
library(mixtools)

files=list.files(path = directory, pattern = "dist", all.files = FALSE,
                 full.names = FALSE, recursive = FALSE,
                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

for (file in files){
  crdsize=as.data.frame(data.table::fread(paste0(directory,file), head=TRUE, stringsAsFactors=FALSE))
  crdsize=as.matrix(crdsize)
  crdsize.log = log10(crdsize-2)
  
  mixmdl = normalmixEM(crdsize.log,mu = c(2.17,4.46))
  thres = qnorm(0.95,mean=mixmdl$mu[1],sd=mixmdl$sigma[1])
  print(thres*1000)
  
  pdf(paste0(out_directory,"mixture_model_",file,".pdf"))
  plot(mixmdl,density = TRUE, w = 2,xlab2="log10(Size)",n=40)
  lines(density(crdsize.log), lty=2, lwd=2)
  abline(v=thres,lty=2,lwd=3)
  text(x=thres+1,y=0.3,paste0("Size threshold:\n",round((10^thres)/1000,1)," kb"))
  dev.off()
  
  filter_crds=crdsize[crdsize > thres*1000]
   
  line=paste0('threshold ',thres*1000, '  crd_all  ',length(crdsize),'  crds_after_thresh ',length(filter_crds) , ' file  ' , file)
  write(line, file = paste0(out_directory,"threshold_methylCRD_.txt"),append = TRUE, sep = " ")
  
}

