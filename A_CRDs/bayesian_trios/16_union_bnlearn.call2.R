#!/usr/bin/env Rscript 

# Clean environment ---------------------------------
rm(list=ls())
gc()


# Packages ---------------------------------

library(bnlearn)
library(dplyr)


# Directories and Data ---------------------------------

# directories for mac =================================
DIR_T='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/triplets_data'
DIR_BN='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/BN'

# directories for cluster =================================
DIR_T='/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_data_shared'
DIR_BN='/home/users/a/avalosma/scratch/12_TRIPLETS/BN_union'


# Main ---------------------------------
data_type='hist'
cell_type_ref1='neut'
cell_type_query2='mono'

data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
window='1000000'

for(data_type in data_types){
  for(cell_type_ref1 in cell_types){
    for(cell_type_query2 in cell_types){
      if (cell_type_ref1 != cell_type_query2){
    
        df_total=data.frame()
        
        name_condition=paste0(data_type,'_',cell_type_query2,'_vs_',cell_type_ref1,'_all' )
        cat(name_condition, '  \n')
        
        all.files <- list.files(path=DIR_T, pattern=paste0(data_type,'_',cell_type_query2,'_vs_',cell_type_ref1 ), 
                                full.names=TRUE, recursive=FALSE)
        
        for (file in all.files){
          
          filename=unlist(strsplit(file, "/"))[9] # 8 for mac, 9 for cluster
          nbr=unlist(strsplit(filename, "_"))[5]
          var=unlist(strsplit(filename, "_"))[6]
          gene=unlist(strsplit(filename, "_"))[7]
          crd=paste0(unlist(strsplit(filename, "_"))[8],
                     '_',unlist(strsplit(filename, "_"))[9],
                     '_',unlist(strsplit(filename, "_"))[10])     
          AF=as.numeric(unlist(strsplit(unlist(strsplit(filename, "_"))[12], ".t"))[1])
          dist=unlist(strsplit(filename, "_"))[11] 
     
          d = read.table(file, stringsAsFactor=FALSE, head=TRUE)
          d <- d[,-1] # remove sample names
          colnames(d)=c("G","V","C")       # V = variant # C = CRD activity  # G = Gene
          
          if(!(is.nan((d$V - mean(d$V))/sd(d$V)))[1]){
            d$V = (d$V - mean(d$V))/sd(d$V) # standardization of V (mean=0, sted=1) 
            
            #list of 4 models to test
            m=rep("", 3)
            m[1]="[V][C|V][G|C]"	#V -> C -> G # Scenario where the variants affects the CRD that affects the G (causal)
            m[2]="[V][G|V][C|G]"	#V -> G -> C # Scenario where the variant affects the Gene and the Gene affects the CRD (reactive)
            m[3]="[V][C|V][G|V]"	#V -> G, V-> C # Scenario where the variants affects independently the CRD and Gene (independent)
            
            #calcuate P(D|G) where G=1,2,3 ## networks score of a particualr graph for a particular data set
            loglik=rep(0, 3)
            for (i in 1:3) {
              net = model2network(m[i])
              loglik[i]=score(net, d, type="bge")
            }
            
            #assume constant prior over the 4 network configurations, i.e. P(G) = 0.25, and compute posterior
            prior=rep(1/3, 3)
            posteriors = exp(loglik - max(loglik)) * prior # scores are in log scale. substraction means division!!!  exp exponential value; exp(x), e to the power of x; x^5 = 2.7^5, opposite of log
            posteriors = posteriors / sum(posteriors) # Division by the sum so that all of the probabilities add to 1!!
            
           #write output
      #      cat ("M1=", signif(posteriors[1],3), " M2=", signif(posteriors[2],3), "M3=", signif(posteriors[3],3), "\n")
            df <- data.frame("name" = name_condition,
                             "nbr" = nbr,
                             "var" = var,
                             "gene" = gene,
                             "crd" = crd,
                             "dist" = dist,
                             "AF" = AF,
                             "M1" = signif(posteriors[1],3),
                             "M2" = signif(posteriors[2],3),
                             "M3" = signif(posteriors[3],3),
                             "L1" = loglik[1],
                             "L2" = loglik[2],
                             "L3" = loglik[3])
            df_total=bind_rows(df_total,df )
          }
        }  
        output=paste0(DIR_BN,'/',name_condition,'.txt')
        write.table(df_total, output, quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")
      }
    }
  }
}
      



