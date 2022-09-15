#!/usr/bin/env Rscript 

library(bnlearn)
library(GenABEL)
#read input file
args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]


## d = read.table(input, stringsAsFactor=FALSE, head=TRUE)
# d=read.csv("test.csv", sep=",", dec=".", header=TRUE,stringsAsFactors=FALSE)
# head(d)
dim(d)

d <- d[,-1]

## col <- names(d)
## triplet <- paste(col,collapse=";")

colnames(d)=c("G","V","C")
d$V = (d$V - mean(d$V))/sd(d$V)

#In statistics, @@standardization@@ is the process of putting different variables on the same scale. 
#This process allows you to compare scores between different types of variables. Typically, to standardize variables, you calculate the mean and standard deviation for a variable. Then, for each observed value of the variable, you subtract the mean and divide by the standard deviation.

# V = variant
# C = CRD activity
# G = Gene

#list of 4 models to test
m=rep("", 3)
## m[1]="[V][T|V][G|T]"	#V -> T -> G # Scenario where the variants affects the TE that affects the G
## m[2]="[V][G|V][T|G]"	#V -> G -> T # Scenario where the variant affects the Gene and the Gene affects the TE
## m[3]="[V][T|V][G|V]"	#V -> G, V-> T # Scenario where the variants affects independently the TE and Gene.

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
cat ("M1=", signif(posteriors[1],3), " M2=", signif(posteriors[2],3), "M3=", signif(posteriors[3],3), "\n")
df <- data.frame("M1" = signif(posteriors[1],3),
				 "M2" = signif(posteriors[2],3),
				 "M3" = signif(posteriors[3],3),
				 "L1" = loglik[1],
				 "L2" = loglik[2],
				 "L3" = loglik[3])

write.table(df, output, quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")


