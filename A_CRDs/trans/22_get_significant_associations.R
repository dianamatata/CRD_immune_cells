################################################################################################################
# GOAL: filter significant trans CRD associations
################################################################################################################

#################################### Packages

library(qvalue)


#################################### Function

get_q_value <- function(trans_data, FDR_threshold, name){
  colnames(trans_data) = c("idx1","chr1","midplace","id1","idx2","chr2","midplace2","id2","corr","pval")
  # Q = qvalue(trans_data$pval, pi0=1)
  Q = qvalue(trans_data$pval)
  trans_data$qval=Q$qvalues
  trans_data_q001 = trans_data[which(trans_data$qval < FDR_threshold), ]
  cat(paste0(name,':  ',length(trans_data_q001[,1]),'  ' ))
  write.table(trans_data_q001, paste0(path_out,'/',name,'.significant_',FDR_threshold,'.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE)
}

get_q_value_2FDRs <- function(trans_data, name){
  colnames(trans_data) = c("idx1","chr1","midplace","id1","idx2","chr2","midplace2","id2","corr","pval")
  # Q = qvalue(trans_data$pval, pi0=1)
  Q = qvalue(trans_data$pval)  
  trans_data$qval=Q$qvalues
  trans_data_q001 = trans_data[which(trans_data$qval < 0.01), ]
  trans_data_q005 = trans_data[which(trans_data$qval < 0.05), ]
  cat(paste0(name,':  ',length(trans_data_q001[,1]),'  ',length(trans_data_q001[,1])/length(trans_data[,1]),'            ' ))
  write.table(trans_data_q001, paste0(path_out,'/',name,'.significant_',0.01,'.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE)
  write.table(trans_data_q005, paste0(path_out,'/',name,'.significant_',0.05,'.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#################################### Folders

path_in = "/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:merged"
path_out = "/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant"


#################################### Main

FDR_threshold=0.05 ## do we want 0.05 too?
file.names <- dir(path_in, pattern =".txt.gz")

for(i in 1:length(file.names)){
  cat(i)
  name= substr(file.names[i],1,nchar(file.names[i])-11)
  trans_data = as.data.frame(data.table::fread(paste0(path_in,'/',file.names[i]), head=FALSE, stringsAsFactors=FALSE))
  cat(name)
  # get_q_value(trans_data, 0.05, name)
  tryCatch( { get_q_value_2FDRs(trans_data, name) }
            , error = paste0('error:',name) )
}


#### ERROR HANDLING: Q value: https://github.com/StoreyLab/qvalue/issues/9
# if smooth.spline error:  be conservative and simply force pi0 = 1: qvalue(p, pi0 = 1) ...
# This is also known as the BH procedure from the function p.adjust.
# for debugging

#################################### Compare values with G
# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0/TRANS% cat test.significant1FDR.txt | wc -l
# result more conservative for us if pi=0 : 301345 (70)g
# otherwise same result
# 70: 159423
# 72: 84691
# 73: 116659
